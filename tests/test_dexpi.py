# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Tests for the exchange of systems with DEXPI 2.0 process models.
"""
import os
import subprocess
import sys
import warnings

import pytest
import biosteam as bst
from numpy.testing import assert_allclose
from biosteam.interop import dexpi
from biosteam.interop.dexpi import model, notes, validation, xmlio

CHEMICALS = ['Water', 'Ethanol', 'AceticAcid', 'EthylAcetate']


def create_system(flowsheet):
    """Create a system with a recycle that covers the mapped unit classes."""
    bst.main_flowsheet.set_flowsheet(flowsheet)
    bst.settings.set_thermo(CHEMICALS, cache=True)
    feed = bst.Stream('feed', Water=80, Ethanol=20, AceticAcid=5, price=0.1)
    P1 = bst.Pump('P1', ins=feed, P=2 * 101325)
    H1 = bst.HXutility('H1', ins=P1-0, T=320)
    R1 = bst.SinglePhaseReactor(
        'R1', ins=H1-0, tau=1., T=320, P=101325,
        reactions=bst.Reaction(
            'Ethanol + AceticAcid -> EthylAcetate + Water', 'AceticAcid', 0.9
        ),
    )
    M1 = bst.Mixer('M1', ins=(R1-0, 'recycle'))
    F1 = bst.Flash('F1', ins=M1-0, outs=('vapor', 'liquid'), T=360, P=101325)
    S1 = bst.Splitter('S1', ins=F1-1, outs=(1-M1, 'product'),
                      split={'Water': 0.2, 'Ethanol': 0.3})
    T1 = bst.StorageTank('T1', ins=S1-1, outs='stored')
    return bst.System.from_units('sys', [P1, H1, R1, M1, F1, S1, T1])


def assert_equivalent(system, copy, rtol=1e-9):
    """Assert that two simulated systems have the same streams."""
    assert [type(i) for i in copy.units] == [type(i) for i in system.units]
    assert [i.ID for i in copy.units] == [i.ID for i in system.units]
    copies = copy.flowsheet.stream
    assert len(copy.streams) == len(system.streams)
    for stream in system.streams:
        other = copies[stream.ID]
        assert other.phases == stream.phases
        assert_allclose(other.mol, stream.mol, rtol=rtol, atol=1e-12)
        assert_allclose(other.T, stream.T, rtol=rtol)
        assert_allclose(other.P, stream.P, rtol=rtol)


def simulate_from_scratch(system):
    """Simulate a system starting from empty process streams, as for a new
    system; converged recycles otherwise depend on the initial guess within
    the convergence tolerance."""
    for stream in system.streams:
        if stream.source:
            stream.empty()
    system.simulate()


@pytest.fixture
def system():
    system = create_system('dexpi_original')
    system.simulate()
    return system


# %% XML layer

def test_xml_round_trip_preserves_all_value_kinds():
    root = model.DexpiObject(model.ENGINEERING_MODEL, data={
        'Text': 'a < b & "c"',
        'Integer': 3,
        'Double': 0.1 + 0.2,
        'Infinity': float('inf'),
        'Boolean': True,
        'Undefined': model.UNDEFINED,
        'Many': ['l', 'g'],
        'Enumeration': model.process_enumeration('PortDirection', 'Inlet'),
        'Quantity': model.physical_quantity(300., 'TemperatureUnit', 'Kelvin'),
    })
    child = root.add('Children', model.DexpiObject(model.process_type('Stream')))
    root.refer('Favorite', child)
    copy = xmlio.loads(xmlio.dumps(model.DexpiDocument(root, name='test'))).root
    assert copy.data['Text'] == 'a < b & "c"'
    assert copy.data['Integer'] == 3 and isinstance(copy.data['Integer'], int)
    assert copy.data['Double'] == 0.1 + 0.2  # Exact, not approximate
    assert copy.data['Infinity'] == float('inf')
    assert copy.data['Boolean'] is True
    assert copy.data['Undefined'] is model.UNDEFINED
    assert copy.get('Undefined', 'default') == 'default'
    assert copy.data['Many'] == ['l', 'g']
    assert copy.data['Enumeration'].literal == 'Inlet'
    assert model.read_quantity(copy.data['Quantity'])[0] == 300.
    assert copy.target('Favorite') is copy.child('Children')


@pytest.mark.parametrize('text, message', [
    ('<Model><Object type="Core/EngineeringModel"', 'not well-formed'),
    ('<PlantModel/>', "expected a 'Model' root element"),
    ('<Model/>', 'exactly one root object'),
    ('<Model><Object id="a" type="Core/EngineeringModel">'
     '<References property="X" objects="#missing"/></Object></Model>',
     'undefined object'),
    ('<!DOCTYPE Model [<!ENTITY x "y">]><Model/>', 'not allowed'),
])
def test_malformed_documents_are_rejected(text, message):
    with pytest.raises(dexpi.DexpiFormatError, match=message):
        xmlio.loads(text)


def test_import_prefixes_are_normalized():
    document = xmlio.loads(
        '<Model name="x" uri="urn:x">'
        '<Import prefix="c" source="https://data.dexpi.org/models/2.0.0/Core.xml"/>'
        '<Import prefix="p" source="https://data.dexpi.org/models/2.0.0/Process.xml"/>'
        '<Object id="a" type="c/EngineeringModel">'
        '<Components property="ConceptualModel">'
        '<Object id="b" type="p/ProcessModel"/></Components></Object></Model>'
    )
    assert document.root.type == model.ENGINEERING_MODEL
    assert document.process_model is not None


# %% Structural validation

def test_validation_reports_dangling_streams_and_ports():
    root = model.DexpiObject(model.ENGINEERING_MODEL)
    process = root.add('ConceptualModel', model.DexpiObject(model.PROCESS_MODEL))
    step = process.add('ProcessSteps', model.DexpiObject(
        model.process_type('Pumping'), data={'Identifier': 'P1'}))
    port = step.add('Ports', model.DexpiObject(
        model.process_type('MaterialPort'),
        data={'Identifier': 'in0',
              'NominalDirection': model.process_enumeration('PortDirection', 'Inlet')}))
    stream = process.add('ProcessConnections', model.DexpiObject(
        model.process_type('Stream'), data={'Identifier': 's1'}))
    stream.refer('Source', port)  # An inlet port cannot be a source
    messages = [i.message for i in validation.find_issues(model.DexpiDocument(root))]
    assert any('is not connected' in i for i in messages)
    assert any('expected an Outlet port' in i for i in messages)
    assert any('has no Target port' in i for i in messages)
    with pytest.raises(dexpi.DexpiValidationError) as error:
        dexpi.from_dexpi(model.DexpiDocument(root))
    assert len(error.value.issues) == len(messages)


def test_document_without_process_model_is_rejected():
    document = model.DexpiDocument(model.DexpiObject(model.ENGINEERING_MODEL))
    with pytest.raises(dexpi.DexpiValidationError, match='ProcessModel'):
        dexpi.from_dexpi(document)
    with pytest.raises(dexpi.DexpiImportError, match='process model'):
        dexpi.from_dexpi(document, validate=False)


# %% Export

def test_export_structure(system):
    document = dexpi.DexpiExporter().export(system)
    assert validation.find_issues(document) == []
    process = document.process_model
    steps = {i.get('Identifier'): i for i in process.children('ProcessSteps')}
    classes = {i: j.class_name for i, j in steps.items()}
    assert classes == {
        'P1': 'Pumping', 'H1': 'ExchangingThermalEnergy',
        'R1': 'ReactingChemicals', 'M1': 'MixingSimple',
        'F1': 'SeparatingByFlash', 'S1': 'SplittingMaterial',
        'T1': 'StoringInTank',
        'feed.feed': 'Source', 'product.vapor': 'Sink', 'product.stored': 'Sink',
    }
    # Every stream is connected on both ends
    streams = process.children('ProcessConnections')
    assert {i.get('Identifier') for i in streams} == {i.ID for i in system.streams}
    assert all(i.target('Source') and i.target('Target') for i in streams)
    # Standard DEXPI attributes
    reactor = steps['R1']
    assert reactor.get('Method').literal == 'Tank'
    assert dexpi.quantities.first(reactor, 'Temperature', 'temperature') == 320.
    assert dexpi.quantities.first(reactor, 'Volume', 'volume') > 0.
    # Components carry database identifiers
    components = process.children('MaterialComponents')
    assert [i.get('Identifier') for i in components] == CHEMICALS
    water = components[0]
    assert water.class_name == 'PureMaterialComponent'
    identifiers = {i.get('Context'): i.get('Value')
                   for i in water.children('PersistentIdentifiers')}
    assert identifiers['CAS'] == '7732-18-5'
    assert identifiers['InChIKey'] == 'XLYOFNOQVPJJNP-UHFFFAOYSA-N'
    # BioSTEAM parameters travel as notes
    note = notes.read_note(steps['S1'])
    assert note['class'] == 'Splitter'
    assert note['parameters']['split']['value'] == {'Water': 0.2, 'Ethanol': 0.3}
    assert notes.read_note(reactor)['parameters']['reactions']['value']['X'] == 0.9


def test_stream_data_uses_dexpi_units(system):
    document = dexpi.DexpiExporter().export(system)
    streams = {i.get('Identifier'): i
               for i in document.process_model.children('ProcessConnections')}
    feed = streams['feed']
    assert dexpi.quantities.first(feed, 'MassFlow', 'mass flow') == pytest.approx(
        system.flowsheet.stream.feed.F_mass)
    state = feed.target('MaterialStateReference').child('State')
    flows = state.child('Composition').children('MoleFlow')
    unit = flows[0].get('Value').data['Unit']
    assert unit.literal == 'KilomolePerSecond'
    assert flows[0].get('Value').data['Value'] == pytest.approx(80. / 3600.)
    assert notes.read_note(feed)['price'] == 0.1


def test_export_to_file_and_method_api(system, tmp_path):
    path = tmp_path / 'system.xml'
    assert system.to_dexpi(path) is None
    assert path.read_text(encoding='utf-8') == system.to_dexpi(
        timestamp=xmlio.read(path).root.get('ExportDateTime'))
    copy = bst.System.from_dexpi(path, flowsheet='dexpi_from_file')
    assert copy.ID == 'sys'
    assert len(copy.units) == len(system.units)


# %% Round trip

def test_round_trip_simulation_parity(system):
    copy = dexpi.from_dexpi(system.to_dexpi(), flowsheet='dexpi_copy')
    assert copy.flowsheet is not system.flowsheet
    assert bst.main_flowsheet.ID == 'dexpi_original'  # Restored after import
    # Streams hold the exported results before simulation...
    assert_equivalent(system, copy, rtol=1e-12)
    # ...and the copy reproduces them when simulated.
    simulate_from_scratch(copy)
    assert_equivalent(system, copy)
    assert_allclose(copy.installed_equipment_cost, system.installed_equipment_cost)
    assert copy.flowsheet.stream.feed.price == 0.1
    assert copy.flowsheet.unit.S1.isplit['Ethanol'] == 0.3


def structure(xml):
    """Return the classes, connectivity and annotations of a process model."""
    process = xmlio.loads(xml).process_model
    steps = {
        i.get('Identifier'): (
            i.class_name,
            sorted(j.get('Identifier') for j in i.children('Ports')),
            notes.read_note(i),
        )
        for i in process.children('ProcessSteps')
    }
    connections = {
        i.get('Identifier'): (
            i.target('Source').get('Identifier'),
            i.target('Target').get('Identifier'),
        )
        for i in process.children('ProcessConnections')
    }
    return steps, connections


def test_second_round_trip_is_stable(system):
    first = system.to_dexpi()
    copy = dexpi.from_dexpi(first, flowsheet='dexpi_copy_1')
    second = copy.to_dexpi()
    assert structure(second) == structure(first)
    again = dexpi.from_dexpi(second, flowsheet='dexpi_copy_2')
    simulate_from_scratch(again)
    assert_equivalent(system, again)


def test_user_chemicals_are_used_as_is(system):
    thermo = bst.settings.get_thermo()
    copy = dexpi.from_dexpi(system.to_dexpi(), flowsheet='dexpi_same_thermo')
    assert copy.units[0].thermo is thermo  # Active package matches the document
    other = bst.Thermo(bst.Chemicals(CHEMICALS, cache=True))
    copy = dexpi.from_dexpi(system.to_dexpi(), flowsheet='dexpi_given_thermo',
                            chemicals=other)
    assert copy.units[0].thermo is other


def test_multiphase_streams_round_trip():
    bst.main_flowsheet.set_flowsheet('dexpi_multiphase')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.MultiStream('feed', l=[('Water', 60), ('Ethanol', 10)],
                           g=[('Water', 5), ('Ethanol', 25)], T=355)
    M1 = bst.Mixer('M1', ins=feed, outs='mixed')
    system = bst.System.from_units('multiphase_sys', [M1])
    system.simulate()
    copy = dexpi.from_dexpi(system.to_dexpi(), flowsheet='dexpi_multiphase_copy')
    feed_copy = copy.flowsheet.stream.feed
    assert set(feed_copy.phases) == {'l', 'g'}
    assert_allclose(feed_copy.imol['g', 'Ethanol'], 25.)
    assert_allclose(feed_copy.imol['l', 'Water'], 60.)


# %% Mapping registry and fallbacks

def test_registry_lookup_follows_inheritance():
    registry = dexpi.default_registry
    class SpecialPump(bst.Pump):
        pass

    assert registry.lookup_unit_class(SpecialPump).dexpi_class == 'Pumping'
    assert registry.lookup_unit_class(bst.MixTank).dexpi_class == 'StoringInTank'
    assert registry.lookup_unit_class(bst.Unit).dexpi_class == 'ProcessStep'
    assert registry.lookup_dexpi_class('RotaryMixing').unit_class is bst.Mixer
    assert registry.lookup_dexpi_class('Cooling').unit_class is bst.HXutility
    assert registry.lookup_dexpi_class('Sieving') is None
    assert registry.resolve_class_name('Pump') is bst.Pump
    assert registry.resolve_class_name('os.system') is None
    assert registry.resolve_class_name('Stream') is None  # Not a unit


def test_custom_registration_does_not_change_default_registry():
    class Centrifuge(bst.Splitter):
        pass

    registry = dexpi.default_registry.copy()
    registry.register(Centrifuge, 'SeparatingByCentrifugalForce', import_default=True)
    assert registry.lookup_dexpi_class('SeparatingByCentrifugalForce').unit_class is Centrifuge
    assert dexpi.default_registry.lookup_dexpi_class('SeparatingByCentrifugalForce') is None
    with pytest.raises(TypeError):
        registry.register(object, 'ProcessStep')


class Booster(bst.Unit):
    """A unit class that is not part of BioSTEAM or the registry."""
    _N_ins = _N_outs = 1

    def _init(self, factor=2.):
        self.factor = factor

    def _run(self):
        self.outs[0].copy_like(self.ins[0])
        self.outs[0].P *= self.factor


def test_unmapped_unit_class_round_trips_through_generic_step(monkeypatch):
    bst.main_flowsheet.set_flowsheet('dexpi_unmapped')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=10, Ethanol=1)
    B1 = Booster('B1', ins=feed, outs='boosted', factor=3.)
    M1 = bst.Mixer('M1', ins=(B1-0, bst.Stream('makeup', Water=1)), outs='mixed')
    system = bst.System.from_units('unmapped_sys', [B1, M1])
    system.simulate()
    xml = system.to_dexpi()
    step = {i.get('Identifier'): i for i in
            xmlio.loads(xml).process_model.children('ProcessSteps')}['B1']
    assert step.class_name == 'ProcessStep'
    assert notes.read_note(step)['parameters'] == {'factor': 3.}

    # Unit classes defined by the user are restored by name
    restored = dexpi.from_dexpi(xml, flowsheet='dexpi_unmapped_restored')
    assert type(restored.units[0]) is Booster
    assert restored.units[0].factor == 3.

    # Where the class is not defined, the topology is still complete
    monkeypatch.delattr(bst, 'Booster')
    with pytest.warns(dexpi.DexpiWarning, match="'Booster'.*not available"):
        copy = dexpi.from_dexpi(xml, flowsheet='dexpi_unmapped_copy')
    B1_copy, M1_copy = copy.units
    assert isinstance(B1_copy, dexpi.GenericProcessStep)
    assert B1_copy.outs[0].sink is M1_copy
    copy.simulate()
    assert_allclose(copy.flowsheet.stream.mixed.mol, system.flowsheet.stream.mixed.mol)

    # A registry can provide the class explicitly
    registry = dexpi.default_registry.copy()
    registry.register(Booster, 'Pumping')
    restored = dexpi.from_dexpi(xml, flowsheet='dexpi_unmapped_registered',
                                registry=registry)
    assert type(restored.units[0]) is Booster


def third_party_document():
    """Return a process model as another tool could write it: no BioSTEAM
    notes, non-SI units, and a class BioSTEAM has no model for."""
    path = os.path.join(os.path.dirname(__file__), 'dexpi_third_party.xml')
    with open(path, encoding='utf-8') as file:
        return file.read()


def test_import_of_third_party_document():
    bst.main_flowsheet.set_flowsheet('dexpi_third_party_main')
    system = dexpi.from_dexpi(third_party_document(), flowsheet='dexpi_third_party',
                              resolver=bst.interop.chemicals.default_resolver())
    pump, screen = system.units
    assert system.ID == 'plant'
    # Mapped by DEXPI class, parameterized from standard attributes
    assert type(pump) is bst.Pump and pump.ID == 'P-101'
    assert pump.P == pytest.approx(4e5)
    # No model for sieving: kept as a generic step with its DEXPI class
    assert isinstance(screen, dexpi.GenericProcessStep)
    assert screen.dexpi_class == 'Sieving'
    # Components are resolved from CAS number and name
    chemicals = pump.thermo.chemicals
    assert chemicals.IDs == ('7732-18-5', 'EtOH')
    assert chemicals.EtOH.CAS == '64-17-5'
    # Stream data are converted to BioSTEAM units
    raw = pump.ins[0]
    assert raw.T == pytest.approx(313.15) and raw.P == pytest.approx(1.5e5)
    assert_allclose(raw.mol, [36., 7.2])
    system.simulate()
    assert_allclose(screen.outs[0].mol, [36., 7.2])
    assert screen.outs[0].P == pytest.approx(4e5)

    # Exporting again restores the original class and attributes
    steps = {i.get('Identifier'): i for i in
             xmlio.loads(system.to_dexpi()).process_model.children('ProcessSteps')}
    assert steps['S-102'].class_name == 'Sieving'
    assert steps['P-101'].class_name == 'Pumping'


def test_step_that_cannot_be_modeled_falls_back_to_generic_step():
    # A heat exchanger without outlet temperature cannot be an HXutility
    text = third_party_document().replace(
        'Process/Process.Pumping', 'Process/Process.Cooling')
    bst.main_flowsheet.set_flowsheet('dexpi_fallback_main')
    with pytest.warns(dexpi.DexpiWarning, match='generic process step'):
        system = dexpi.from_dexpi(text, flowsheet='dexpi_fallback',
                                  resolver=bst.interop.chemicals.default_resolver())
    assert [type(i) for i in system.units] == [dexpi.GenericProcessStep] * 2
    assert system.units[0].dexpi_class == 'Cooling'


def test_nonserializable_parameters_are_reported(system):
    reactor = system.flowsheet.unit.R1
    reactor.reactions = bst.Reaction(
        'Ethanol + AceticAcid -> EthylAcetate + Water', 'AceticAcid',
        rate=lambda stream: 0.1,
    )
    with pytest.warns(bst.interop.reactions.ReactionSerializationWarning):
        xml = system.to_dexpi()
    with pytest.warns(bst.interop.reactions.ReactionSerializationWarning):
        copy = dexpi.from_dexpi(xml, flowsheet='dexpi_no_kinetics')
    assert copy.flowsheet.unit.R1.reactions.X == 0.


# %% Optional dependencies

def test_biosteam_import_does_not_load_interop_or_optional_packages():
    code = (
        "import sys, biosteam\n"
        "loaded = [i for i in sys.modules if i.split('.')[0] in "
        "('pydexpi', 'pyenzyme') or i.startswith('biosteam.interop')]\n"
        "assert not loaded, loaded\n"
        "from biosteam.interop import dexpi, enzymeml\n"
        "loaded = [i for i in sys.modules if i.split('.')[0] in ('pydexpi', 'pyenzyme')]\n"
        "assert not loaded, loaded\n"
    )
    result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


# %% Conformance with an independent DEXPI implementation

def test_export_loads_in_pydexpi_strict_mode(system):
    pydexpi = pytest.importorskip('pydexpi')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        loaded = pydexpi.load_xml_string(system.to_dexpi(), strict=True)
    process = loaded.conceptualModel
    assert type(process).__name__ == 'ProcessModel'
    steps = {i.identifier: i for i in process.processSteps}
    assert type(steps['R1']).__name__ == 'ReactingChemicals'
    assert type(steps['P1']).__name__ == 'Pumping'
    assert len(process.processConnections) == len(system.streams)
    feed = next(i for i in process.processConnections if i.identifier == 'feed')
    assert feed.temperature[0].value.value == pytest.approx(298.15)
    assert len(feed.materialTemplateReference.listOfComponents.component) == 4


def test_document_written_by_pydexpi_is_imported(tmp_path):
    pytest.importorskip('pydexpi')
    from pydexpi.loaders.xml_serializer import DexpiXmlSerializer
    from pydexpi.pydexpi_classes import pydexpi_classes as d

    def quantity(value, unit):
        return [d.QualifiedValue(value=d.PhysicalQuantity(value=value, unit=unit))]

    water = d.PureMaterialComponent(identifier='Water', label='water')
    components = d.ListOfMaterialComponents(component=[water])
    text = d.MultiLanguageString(
        singleLanguageStrings=[d.SingleLanguageString(value='water')])
    template = d.MaterialTemplate(
        listOfComponents=components, description=text, identifier='T',
        label='T', numberOfMaterialComponents=1, numberOfPhases=1)
    composition = d.Composition(
        moleFlow=quantity(0.5, d.MoleFlowRateUnit.KILOMOLEPERSECOND))
    state = d.MaterialState(
        state=d.MaterialStateType(composition=composition, description=text,
                                  identifier='total', label='total'),
        description='state', identifier='state', label='state')

    def stream(identifier, **kwargs):
        return d.Stream.model_construct(
            identifier=identifier, materialTemplateReference=template, **kwargs)

    def port(identifier, direction, connection):
        port = d.MaterialPort.model_construct(
            identifier=identifier, nominalDirection=direction,
            connectorReference=connection)
        if direction is d.PortDirection.OUTLET:
            connection.source = port
        else:
            connection.target = port
        return port

    inlet = stream('inlet', materialStateReference=state,
                   temperature=quantity(25., d.TemperatureUnit.DEGREECELSIUS))
    outlet = stream('outlet')
    IN, OUT = d.PortDirection.INLET, d.PortDirection.OUTLET
    steps = [
        d.Source(identifier='source', ports=[port('out0', OUT, inlet)]),
        d.Pumping(identifier='P1', ports=[port('in0', IN, inlet), port('out0', OUT, outlet)],
                  pressure=quantity(3., d.PressureAbsoluteUnit.BAR)),
        d.Sink(identifier='sink', ports=[port('in0', IN, outlet)]),
    ]
    process = d.ProcessModel(
        processSteps=steps, processConnections=[inlet, outlet],
        materialComponents=[water], listsOfMaterialComponents=[components],
        materialTemplates=[template], materialStates=[state])
    DexpiXmlSerializer().save(d.EngineeringModel(conceptualModel=process),
                              tmp_path, 'pydexpi_model')

    bst.main_flowsheet.set_flowsheet('dexpi_from_pydexpi_main')
    system = dexpi.from_dexpi(tmp_path / 'pydexpi_model.xml',
                              flowsheet='dexpi_from_pydexpi',
                              resolver=bst.interop.chemicals.default_resolver())
    pump, = system.units
    assert type(pump) is bst.Pump and pump.P == pytest.approx(3e5)
    assert pump.ins[0].imol['Water'] == pytest.approx(1800.)
    assert pump.ins[0].T == pytest.approx(298.15)
    system.simulate()
    assert pump.outs[0].P == pytest.approx(3e5)
