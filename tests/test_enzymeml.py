# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Tests for loading enzyme kinetics from EnzymeML documents.
"""
import builtins
import json
import math
import warnings

import pytest
import biosteam as bst
from numpy.testing import assert_allclose
from biosteam.interop import enzymeml, examples
from biosteam.interop.enzymeml import RateLaw, RateLawError, SteadyStateCSTRRate


def unit(*base_units):
    return {'base_units': [
        {'kind': kind, 'exponent': exponent, 'multiplier': multiplier, 'scale': scale}
        for kind, exponent, multiplier, scale in base_units
    ]}


mM = unit(('mole', 1, 1., -3), ('litre', -1, 1., 0))
per_second = unit(('second', -1, 1., 0))


def michaelis_menten_document(**reaction):
    """An EnzymeML document of S -> P with Michaelis-Menten kinetics."""
    return {
        'name': 'Michaelis-Menten',
        'small_molecules': [
            {'id': 's0', 'name': 'ethanol',
             'references': ['https://commonchemistry.cas.org/detail?cas_rn=64-17-5']},
            {'id': 's1', 'name': 'acetaldehyde', 'inchikey': 'IKHGUXGNUITLKF-UHFFFAOYSA-N'},
        ],
        'proteins': [{'id': 'p0', 'name': 'enzyme', 'sequence': 'GAG'}],
        'reactions': [{
            'id': 'r0', 'name': 'oxidation',
            'reactants': [{'species_id': 's0', 'stoichiometry': 1}],
            'products': [{'species_id': 's1', 'stoichiometry': 1}],
            'modifiers': [{'species_id': 'p0', 'role': 'biocatalyst'}],
            'kinetic_law': {'equation': 'k_cat * p0 * s0 / (K_M + s0)',
                            'equation_type': 'rateLaw'},
            **reaction,
        }],
        'parameters': [
            {'id': 'k_cat', 'symbol': 'k_cat', 'value': 10., 'unit': per_second},
            {'id': 'K_M', 'symbol': 'K_M', 'value': 2., 'unit': mM},
        ],
        'measurements': [{'id': 'm0', 'ph': 7., 'temperature': 25.,
                          'temperature_unit': unit(('celsius', 1, 1., 0))}],
    }


# %% Reading

def test_read_reference_document():
    document = enzymeml.read(examples.FDH_ENZYMEML)
    assert set(document.species) == {'formate', 'nad', 'co2', 'nadh', 'fdh'}
    assert document.species['fdh'].kind == 'protein'
    assert document.species['fdh'].ecnumber == '1.17.1.9'
    reaction = document.get_reaction()
    assert reaction.stoichiometry == {'formate': -1., 'nad': -1., 'co2': 1., 'nadh': 1.}
    assert reaction.modifiers == {'fdh': 'biocatalyst'}
    assert document.temperature == 303.15 and document.pH == 7.5
    # Parameters are converted to kmol, m3 and hr
    assert_allclose(document.parameters['kcat'].coherent_value, 3.7 * 3600)
    assert_allclose(document.parameters['K_m_formate'].coherent_value, 5.6e-3)


def test_read_accepts_dictionaries_paths_and_pydantic_like_objects(tmp_path):
    data = michaelis_menten_document()
    path = tmp_path / 'document.json'
    path.write_text(json.dumps(data))

    class PydanticLike:
        def model_dump(self, mode=None):
            return data

    documents = [enzymeml.read(i) for i in (data, path, str(path), PydanticLike())]
    assert all(i == documents[0] for i in documents)
    assert enzymeml.read(documents[0]) is documents[0]
    assert documents[0].temperature == pytest.approx(298.15)


def test_rate_law_from_equations_list():
    data = michaelis_menten_document(kinetic_law=None)
    data['equations'] = [{'species_id': 'r0', 'equation': 'k_cat * p0',
                          'equation_type': 'rateLaw'}]
    assert enzymeml.read(data).reactions['r0'].rate_law == 'k_cat * p0'


@pytest.mark.parametrize('change, message', [
    (lambda d: d['reactions'][0]['reactants'].append({'species_id': 'x'}),
     "undefined species 'x'"),
    (lambda d: d['small_molecules'].append({'id': 's0'}), 'duplicate species'),
    (lambda d: d['small_molecules'].append({'name': 'no id'}), "without 'id'"),
    (lambda d: d['parameters'][0].update(value='fast'), 'expected a number'),
])
def test_invalid_documents_are_rejected(change, message):
    data = michaelis_menten_document()
    change(data)
    with pytest.raises(enzymeml.EnzymeMLError, match=message):
        enzymeml.read(data)


def test_invalid_json_is_rejected(tmp_path):
    path = tmp_path / 'broken.json'
    path.write_text('{not json')
    with pytest.raises(enzymeml.EnzymeMLError, match='not an EnzymeML JSON document'):
        enzymeml.read(path)
    with pytest.raises(enzymeml.EnzymeMLError, match='JSON object'):
        enzymeml.EnzymeMLDocument.from_dict([1, 2])


def test_omex_without_pyenzyme_raises_helpful_error(monkeypatch, tmp_path):
    real_import = builtins.__import__

    def no_pyenzyme(name, *args, **kwargs):
        if name.split('.')[0] == 'pyenzyme':
            raise ImportError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', no_pyenzyme)
    with pytest.raises(ImportError, match=r'biosteam\[enzymeml\]'):
        enzymeml.read(tmp_path / 'experiment.omex')


def test_reference_document_is_valid_enzymeml():
    pyenzyme = pytest.importorskip('pyenzyme')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        document = pyenzyme.read_enzymeml(examples.FDH_ENZYMEML)
    assert len(document.reactions) == 1
    # A pyenzyme document is accepted directly and reads the same
    assert enzymeml.read(document).reactions == enzymeml.read(examples.FDH_ENZYMEML).reactions
    assert enzymeml.read(document).parameters == enzymeml.read(examples.FDH_ENZYMEML).parameters


# %% Units

def test_unit_conversion():
    factor = enzymeml.conversion_factor
    assert factor(None) == 1.
    assert_allclose(factor(mM), 1e-3)  # mmol/L = mol/m3 = 1e-3 kmol/m3
    assert_allclose(factor(per_second), 3600.)
    assert_allclose(factor(unit(('mole', 1, 1., -6), ('litre', -1, 1., -3))), 1e-3)  # umol/mL
    assert_allclose(factor(unit(('second', -1, 60., 0))), 60.)  # 1/min
    # Minutes as encoded by pyenzyme (multiplier 60 with scale 1)
    assert_allclose(factor(unit(('second', -1, 60., 1))), 60.)
    assert_allclose(factor(unit(('mole', 1, 1., -3), ('litre', -1, 1., 0), ('second', -1, 1., 0))), 3.6)
    with pytest.raises(enzymeml.UnitError, match="'volt'"):
        factor(unit(('volt', 1, 1., 0)))


# %% Rate laws

def test_rate_law_evaluation():
    law = RateLaw('vmax * S^2 / (K**2 + S^2) + exp(-E / 2) - min(1, 2)')
    assert law.symbols == {'vmax', 'S', 'K', 'E'}
    expected = 3. * 4. / (1. + 4.) + math.exp(-1.) - 1.
    assert law({'vmax': 3., 'S': 2., 'K': 1., 'E': 2.}) == pytest.approx(expected)


@pytest.mark.parametrize('equation', [
    '__import__("os").system("true")',
    'S.__class__',
    'S[0]',
    'lambda: 1',
    '[i for i in S]',
    'S if S else 1',
    '"text" * 3',
    'open("file")',
    'exp(x=1)',
    'S; S',
    '',
])
def test_unsafe_or_malformed_rate_laws_are_rejected(equation):
    with pytest.raises(RateLawError):
        RateLaw(equation)


def test_rate_law_errors_are_explicit():
    with pytest.raises(RateLawError, match="undefined names \\['K'\\]"):
        RateLaw('S / K', names={'S'})
    law = RateLaw('S / K')
    with pytest.raises(RateLawError, match="no value for 'K'"):
        law({'S': 1.})
    with pytest.raises(RateLawError, match='could not be evaluated'):
        law({'S': 1., 'K': 0.})
    with pytest.raises(RateLawError, match='could not be evaluated'):
        RateLaw('9 ** 9 ** 9 ** 9')({})  # Overflows instead of hanging


# %% Chemicals

def test_species_identities():
    document = enzymeml.read(michaelis_menten_document())
    ethanol, acetaldehyde, enzyme = document.chemical_identities()
    assert ethanol.CAS == '64-17-5' and ethanol.pure
    assert acetaldehyde.InChIKey == 'IKHGUXGNUITLKF-UHFFFAOYSA-N'
    assert not enzyme.pure and enzyme.phase == 'l'
    assert enzyme.MW == pytest.approx(203.2, abs=0.01)  # From the sequence
    fdh = enzymeml.read(examples.FDH_ENZYMEML).species['nadh'].to_identity('NADH')
    assert fdh.ID == 'NADH' and fdh.ChEBI == 'CHEBI:16908'


def test_create_chemicals_with_renamed_species():
    document = enzymeml.read(michaelis_menten_document()).with_chemical_IDs(
        {'s0': 'Ethanol', 's1': 'Acetaldehyde', 'p0': 'Enzyme'})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        chemicals = document.create_chemicals()
    assert chemicals.IDs == ('Water', 'Ethanol', 'Acetaldehyde', 'Enzyme')
    assert chemicals.Ethanol.CAS == '64-17-5'
    assert chemicals.Acetaldehyde.CAS == '75-07-0'
    assert chemicals.Enzyme.MW == pytest.approx(203.2, abs=0.01)
    with pytest.raises(enzymeml.EnzymeMLError, match="no species \\['x'\\]"):
        document.with_chemical_IDs({'x': 'X'})


# %% Kinetics

@pytest.fixture
def michaelis_menten():
    bst.main_flowsheet.set_flowsheet('enzymeml_michaelis_menten')
    document = enzymeml.read(michaelis_menten_document()).with_chemical_IDs(
        {'s0': 'Ethanol', 's1': 'Acetaldehyde', 'p0': 'Enzyme'})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        bst.settings.set_thermo(document.create_chemicals())
    return document


def test_cstr_conversion_matches_analytical_solution(michaelis_menten):
    document = michaelis_menten
    tau = 0.5  # hr
    feed = bst.Stream('feed', Water=1000., units='kg/hr')
    feed.imol['Ethanol'] = 0.01
    feed.imol['Enzyme'] = 1e-7
    reactor = bst.SinglePhaseReactor('R1', ins=feed, tau=tau, T=298.15, P=101325)
    reaction = document.parameterize(reactor)
    assert reactor.reactions is reaction
    reactor.run()

    # Steady-state balance of a CSTR with Michaelis-Menten kinetics:
    # (S0 - S) / tau = vmax * S / (K_M + S), a quadratic equation in S.
    flow = feed.F_vol  # m3/hr
    S0 = 0.01 / flow  # kmol/m3
    vmax = 10. * 3600. * (1e-7 / flow)  # kmol/m3/hr
    K_M = 2e-3  # kmol/m3
    b = K_M + vmax * tau - S0
    S = (-b + math.sqrt(b * b + 4. * K_M * S0)) / 2.
    effluent = reactor.outs[0]
    assert_allclose(effluent.imol['Ethanol'] / flow, S, rtol=1e-10)
    assert_allclose(effluent.imol['Acetaldehyde'], 0.01 - effluent.imol['Ethanol'], rtol=1e-10)
    assert effluent.imol['Enzyme'] == 1e-7  # The catalyst is conserved

    # The kinetics follow the residence time of the reactor they are bound to
    conversion = 1. - effluent.imol['Ethanol'] / 0.01
    reactor.tau = 2.
    reactor.run()
    assert 1. - reactor.outs[0].imol['Ethanol'] / 0.01 > conversion


def test_cstr_rate_limits(michaelis_menten):
    reaction = michaelis_menten.to_reaction(residence_time=1.)
    rate = reaction.rate
    stream = bst.Stream(None, Water=1000., units='kg/hr')
    assert rate(stream) == 0.  # No substrate
    stream.imol['Ethanol'] = 0.01
    assert rate(stream) == 0.  # No enzyme
    stream.imol['Enzyme'] = 1.  # Vast excess of enzyme: complete conversion
    assert rate(stream) == pytest.approx(0.01, rel=1e-6)
    rate.residence_time = None
    with pytest.raises(RuntimeError, match='residence time is not defined'):
        rate(stream)
    with pytest.raises(TypeError, match="no residence time attribute 'tau'"):
        rate.bind(object())


def test_reversible_rate_law_reaches_equilibrium():
    bst.main_flowsheet.set_flowsheet('enzymeml_reversible')
    bst.settings.set_thermo(['Water', 'Ethanol', 'Propanol'], cache=True)
    rate = SteadyStateCSTRRate(
        'k * (A - B / K_eq)', {'k': 1e6, 'K_eq': 3.},
        species={'A': 'Ethanol', 'B': 'Propanol'},
        stoichiometry={'Ethanol': -1., 'Propanol': 1.},
        reactant='Ethanol', residence_time=1.,
    )
    forward = bst.Stream(None, Water=100., Ethanol=4.)
    assert rate(forward) == pytest.approx(3., rel=1e-5)  # B / A = 3 at equilibrium
    reverse = bst.Stream(None, Water=100., Propanol=4.)
    assert rate(reverse) == pytest.approx(-1., rel=1e-5)


def test_stoichiometric_numbers_scale_the_rate():
    bst.main_flowsheet.set_flowsheet('enzymeml_stoichiometry')
    bst.settings.set_thermo(['Water', 'Ethanol', 'Propanol'], cache=True)
    # 2 A -> B with r = k * A (per reaction event): A is consumed at 2 r
    rate = SteadyStateCSTRRate(
        'k * A', {'k': 1.}, species={'A': 'Ethanol'},
        stoichiometry={'Ethanol': -2., 'Propanol': 1.},
        reactant='Ethanol', residence_time=0.5,
    )
    stream = bst.Stream(None, Water=100., Ethanol=1.)
    # A0 - A = 2 k tau A  =>  A = A0 / 2
    assert rate(stream) == pytest.approx(0.5)
    reaction = bst.Reaction({'Ethanol': -2., 'Propanol': 1.}, 'Ethanol', rate=rate)
    reaction(stream)
    assert_allclose([stream.imol['Ethanol'], stream.imol['Propanol']], [0.5, 0.25])
    with pytest.raises(ValueError, match='not a reactant'):
        SteadyStateCSTRRate('k', {'k': 1.}, {}, {'Ethanol': -1., 'Propanol': 1.},
                            reactant='Propanol')


def test_reaction_requires_kinetic_information(michaelis_menten):
    data = michaelis_menten_document(kinetic_law=None)
    with pytest.raises(enzymeml.EnzymeMLError, match='has no rate law'):
        enzymeml.read(data).to_reaction()
    data = michaelis_menten_document()
    data['parameters'][0]['value'] = None
    with pytest.raises(enzymeml.EnzymeMLError, match="'k_cat'.*has no value"):
        enzymeml.read(data).to_reaction()
    data = michaelis_menten_document()
    data['reactions'][0]['kinetic_law']['equation'] = 'k_cat * unknown'
    with pytest.raises(RateLawError, match='undefined names'):
        enzymeml.read(data).to_reaction()
    with pytest.raises(enzymeml.EnzymeMLError, match="no reaction 'r9'"):
        michaelis_menten.to_reaction('r9')
    with pytest.raises(enzymeml.EnzymeMLError, match='not a reactant'):
        michaelis_menten.to_reaction(reactant='s1')


def test_rate_serialization_round_trip(michaelis_menten):
    reaction = michaelis_menten.to_reaction(residence_time=0.5)
    data = bst.interop.reactions.reaction_to_dict(reaction)
    assert json.loads(json.dumps(data)) == data  # Plain JSON
    rate = data['rate']
    assert rate['type'] == 'SteadyStateCSTRRate'
    assert rate['parameters'] == {'K_M': pytest.approx(2e-3), 'k_cat': pytest.approx(36000.)}
    assert rate['species'] == {'p0': 'Enzyme', 's0': 'Ethanol'}
    assert rate['residence_time'] == 0.5
    assert rate['source']['reaction'] == 'r0'
    assert rate['source']['enzymes'] == [{'id': 'p0', 'name': 'enzyme'}]
    copy = bst.interop.reactions.reaction_from_dict(data, bst.settings.chemicals)
    stream = bst.Stream(None, Water=1000., units='kg/hr')
    stream.imol['Ethanol'] = 0.01
    stream.imol['Enzyme'] = 1e-7
    assert copy.rate(stream) == reaction.rate(stream) > 0.
