# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Integration test of the interoperability workflow: an EnzymeML reaction
model is simulated in a CSTR, exported as a DEXPI process model, imported
again and simulated with the same results.

The reference case is the regeneration of NADH from NAD+ by formate
dehydrogenase (FDH), with formate as sacrificial cosubstrate.
"""
import json
import subprocess
import sys
import warnings

import pytest
import biosteam as bst
from numpy.testing import assert_allclose
from biosteam.interop import chemicals, dexpi, enzymeml, examples
from biosteam.interop.dexpi import notes, xmlio

IDs = ('Water', 'formate', 'nad', 'co2', 'nadh', 'fdh')


@pytest.fixture
def fdh_system():
    bst.main_flowsheet.set_flowsheet('fdh_original')
    bst.settings.set_thermo(examples.create_fdh_chemicals())
    system = examples.create_fdh_system()
    system.simulate()
    return system


def results(system):
    """Return the results of a simulated system as plain data."""
    return {
        'streams': {
            i.ID: {'T': i.T, 'P': i.P, 'phase': i.phase,
                   'mol': [float(j) for j in i.imol[IDs]]}
            for i in system.streams
        },
        'reactor volume': (
            system.flowsheet.unit.R1.design_results['Reactor volume']),
        'duty': system.flowsheet.unit.R1.net_duty,
        'installed cost': system.installed_equipment_cost,
    }


def assert_same_results(actual, expected, rtol=1e-9):
    assert actual['streams'].keys() == expected['streams'].keys()
    for ID, stream in expected['streams'].items():
        other = actual['streams'][ID]
        assert other['phase'] == stream['phase']
        assert_allclose(other['mol'], stream['mol'], rtol=rtol, atol=1e-15)
        assert_allclose([other['T'], other['P']], [stream['T'], stream['P']],
                        rtol=rtol)
    for name in ('reactor volume', 'duty', 'installed cost'):
        assert_allclose(actual[name], expected[name], rtol=rtol)


def test_chemicals_are_resolved_from_enzymeml_identifiers():
    resolved = examples.create_fdh_chemicals()
    assert resolved.IDs == IDs
    assert resolved.formate.CAS == '64-18-6'  # By InChIKey
    assert resolved.nad.CAS == '53-84-9'
    assert resolved.co2.CAS == '124-38-9'
    assert resolved.nadh.formula == 'C21H29N7O14P2'  # From the property record
    assert resolved.fdh.MW == 40400.
    assert 'CHEBI:16908' in resolved.nadh.aliases
    # Without property records, unknown species are placeholders
    document = enzymeml.read(examples.FDH_ENZYMEML)
    with pytest.warns(chemicals.PlaceholderWarning,
                      match="'nadh' was not found"):
        document.create_chemicals()
    with pytest.raises(chemicals.ChemicalResolutionError, match="'nadh'"):
        document.create_chemicals(
            chemicals.default_resolver(placeholders=False))


def test_fdh_cstr_simulation(fdh_system):
    stream = fdh_system.flowsheet.stream
    feed, product = stream.feed, stream.product
    reactor = fdh_system.flowsheet.unit.R1
    converted = feed.imol['nad'] - product.imol['nad']
    assert 0.85 < converted / feed.imol['nad'] < 0.95

    # Stoichiometry: formate + NAD+ -> CO2 + NADH, enzyme conserved
    assert_allclose(feed.imol['formate'] - product.imol['formate'], converted,
                    rtol=1e-9)
    assert_allclose(product.imol['nadh'], converted, rtol=1e-9)
    assert_allclose(product.imol['co2'], converted, rtol=1e-9)
    assert product.imol['fdh'] == feed.imol['fdh']
    # Balanced in atoms
    assert_allclose(product.F_mass, feed.F_mass, rtol=1e-9)
    assert product.T == pytest.approx(303.15) and reactor.T == 303.15

    # The outlet satisfies the CSTR design equation with the rate law and
    # parameters of the EnzymeML document (ordered bi-bi, units of mM and s).
    flow = reactor.ins[0].F_vol  # m3/hr
    nad, formate, fdh = (
        product.imol[i] / flow * 1e3 for i in ('nad', 'formate', 'fdh')
    )
    kcat, K_m_formate, K_m_nad, K_i_nad = 3.7, 5.6, 0.045, 0.04
    rate = kcat * fdh * nad * formate / (
        K_i_nad * K_m_formate + K_m_formate * nad + K_m_nad * formate
        + nad * formate
    )  # mM/s = mol/m3/s
    volume = reactor.tau * flow  # m3
    assert_allclose(rate * 3.6 * volume, converted, rtol=1e-9)  # kmol/hr

    # More enzyme or residence time increases the conversion
    reactor.tau = 4.
    fdh_system.simulate()
    assert feed.imol['nad'] - product.imol['nad'] > converted


def test_fdh_dexpi_export(fdh_system):
    document = xmlio.loads(fdh_system.to_dexpi())
    process = document.process_model
    steps = {i.get('Identifier'): i for i in process.children('ProcessSteps')}
    assert steps['R1'].class_name == 'ReactingChemicals'
    assert steps['P1'].class_name == 'Pumping'
    assert steps['H1'].class_name == 'ExchangingThermalEnergy'
    components = {
        i.get('Identifier'): i for i in process.children('MaterialComponents')
    }
    assert components['formate'].class_name == 'PureMaterialComponent'
    assert components['formate'].get('ChEBI_Identifier') == 'CHEBI:30751'
    assert components['fdh'].class_name == 'CustomMaterialComponent'
    # The kinetics and their provenance are part of the reactor annotation
    rate = notes.read_note(
        steps['R1'])['parameters']['reactions']['value']['rate']
    assert rate['type'] == 'SteadyStateCSTRRate'
    assert rate['parameters']['kcat'] == pytest.approx(3.7 * 3600)
    assert rate['source']['format'] == 'EnzymeML'
    assert rate['source']['enzymes'][0]['ecnumber'] == '1.17.1.9'
    assert 'residence_time' not in rate  # Follows the reactor


def test_fdh_round_trip_with_independent_property_package(fdh_system):
    expected = results(fdh_system)
    xml = fdh_system.to_dexpi()

    # Hydrate the chemicals again from the identifiers in the document,
    # as a user without the original property package would.
    resolver = chemicals.default_resolver(properties=examples.FDH_PROPERTIES)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', chemicals.PlaceholderWarning)
        copy = bst.System.from_dexpi(
            xml, flowsheet='fdh_copy', resolver=resolver)
    assert copy.flowsheet is not fdh_system.flowsheet
    thermo = copy.units[0].thermo
    assert thermo is not bst.settings.get_thermo()
    assert thermo.chemicals.IDs == IDs
    assert [type(i) for i in copy.units] == [type(i) for i in fdh_system.units]

    reactor = copy.flowsheet.unit.R1
    rate = reactor.reactions.rate
    assert isinstance(rate, enzymeml.SteadyStateCSTRRate)
    assert rate.residence_time == reactor.tau == 2.
    assert reactor.reactions.chemicals is thermo.chemicals

    for stream in copy.streams:
        if stream.source:
            stream.empty()
    copy.simulate()
    assert_same_results(results(copy), expected)

    # The copy is a live model: its kinetics follow its own reactor
    reactor.tau = 0.5
    copy.simulate()
    stream = copy.flowsheet.stream
    assert (
        stream.product.imol['nad'] > expected['streams']['product']['mol'][2]
    )


def test_fdh_round_trip_without_any_property_data(fdh_system):
    # Formula, molecular weight and phase travel with the document, so the
    # mass balance is reproduced even if nothing else is known.
    expected = results(fdh_system)
    with pytest.warns(chemicals.PlaceholderWarning,
                      match="'fdh' was not found"):
        copy = dexpi.from_dexpi(fdh_system.to_dexpi(), flowsheet='fdh_bare',
                                resolver=chemicals.default_resolver())
    copy.simulate()
    actual = results(copy)
    for ID, stream in expected['streams'].items():
        assert_allclose(actual['streams'][ID]['mol'], stream['mol'],
                        rtol=1e-9, atol=1e-15)


def test_fdh_kinetics_can_be_replaced_from_enzymeml_after_import(fdh_system):
    copy = dexpi.from_dexpi(fdh_system.to_dexpi(),
                            flowsheet='fdh_reparameterized')
    data = json.loads(open(examples.FDH_ENZYMEML, encoding='utf-8').read())
    for parameter in data['parameters']:
        if parameter['id'] == 'kcat':
            parameter['value'] *= 0.5  # A less active enzyme variant
    reactor = copy.flowsheet.unit.R1
    enzymeml.read(data).parameterize(reactor)
    copy.simulate()
    stream = copy.flowsheet.stream
    slower = 1 - stream.product.imol['nad'] / stream.feed.imol['nad']
    original = fdh_system.flowsheet.stream
    assert slower < (
        1 - original.product.imol['nad'] / original.feed.imol['nad'])


def test_fdh_export_conforms_to_pydexpi(fdh_system):
    pydexpi = pytest.importorskip('pydexpi')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        loaded = pydexpi.load_xml_string(fdh_system.to_dexpi(), strict=True)
    steps = {i.identifier: i for i in loaded.conceptualModel.processSteps}
    reactor = steps['R1']
    assert type(reactor).__name__ == 'ReactingChemicals'
    assert reactor.method.value == 'Tank'
    assert reactor.temperature[0].value.value == pytest.approx(303.15)
    components = {
        i.identifier: i for i in loaded.conceptualModel.materialComponents
    }
    assert components['nad'].chEBI_identifier == 'CHEBI:15846'
    assert type(components['fdh']).__name__ == 'CustomMaterialComponent'


_EXPORT = '''
import json, sys, warnings
warnings.simplefilter('ignore')
import biosteam as bst
from biosteam.interop import examples
sys.path.insert(0, {tests!r})
from test_fdh_cstr_roundtrip import results
bst.settings.set_thermo(examples.create_fdh_chemicals())
system = examples.create_fdh_system()
system.simulate()
system.to_dexpi({xml!r})
json.dump(results(system), open({results!r}, 'w'))
'''

_IMPORT = '''
import json, sys, warnings
warnings.simplefilter('ignore')
import biosteam as bst
sys.path.insert(0, {tests!r})
from test_fdh_cstr_roundtrip import results
system = bst.System.from_dexpi({xml!r})
system.simulate()
json.dump(results(system), open({results!r}, 'w'))
'''


@pytest.mark.slow
def test_fdh_round_trip_between_processes(tmp_path):
    """Export in one Python process and import in another, which shares
    nothing with the first but the DEXPI file."""
    import os
    paths = {
        'tests': os.path.dirname(__file__),
        'xml': str(tmp_path / 'fdh.xml'),
    }
    outputs = []
    for name, script in (('exported', _EXPORT), ('imported', _IMPORT)):
        output = str(tmp_path / f'{name}.json')
        subprocess.run(
            [sys.executable, '-c', script.format(results=output, **paths)],
            check=True, capture_output=True, text=True,
        )
        outputs.append(json.load(open(output)))
    exported, imported = outputs
    assert_same_results(imported, exported)
