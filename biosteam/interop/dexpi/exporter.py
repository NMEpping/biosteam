# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Export of BioSTEAM systems to DEXPI process models.

"""
from __future__ import annotations

from datetime import datetime, timezone
from typing import Optional

import thermosteam as tmo

from ..chemicals import ChemicalIdentity
from . import quantities
from .codecs import default_registry
from .model import (
    ENGINEERING_MODEL,
    PERSISTENT_IDENTIFIER,
    PROCESS_MODEL,
    DexpiDocument,
    DexpiObject,
    multi_language_string,
    process_enumeration,
    process_type,
    qualified_value,
)
from .notes import create_note
from .registry import UnitRegistry

__all__ = ('DexpiExporter',)

_INLET = process_enumeration('PortDirection', 'Inlet')
_OUTLET = process_enumeration('PortDirection', 'Outlet')

#: Identifier contexts of chemicals, by field of a chemical identity.
_IDENTIFIER_CONTEXTS = {
    'CAS': 'CAS',
    'InChIKey': 'InChIKey',
    'InChI': 'InChI',
    'PubChem': 'PubChem',
    'SMILES': 'SMILES',
}


def _quantity(value: float, kind: str) -> DexpiObject:
    return qualified_value(quantities.to_dexpi(value, kind))


class DexpiExporter:
    """
    Translate a BioSTEAM system into a DEXPI document.

    Parameters
    ----------
    registry :
        Mapping between unit classes and DEXPI classes. Defaults to
        :data:`~biosteam.interop.dexpi.default_registry`.
    timestamp :
        Export time written to the document. Defaults to the current time;
        pass a fixed value for reproducible output.

    """

    def __init__(self, registry: Optional[UnitRegistry] = None,
                 timestamp: Optional[datetime] = None):
        self.registry = registry or default_registry
        self.timestamp = timestamp

    def export(self, system) -> DexpiDocument:
        """Return the DEXPI document describing `system`."""
        return _Export(self.registry, self.timestamp, system).document


class _Export:
    """The state of a single export."""

    def __init__(self, registry, timestamp, system):
        import biosteam as bst
        self.registry = registry
        self.system = system
        self.units = list(system.units)
        self._unit_set = set(self.units)
        self.model = model = DexpiObject(PROCESS_MODEL)
        self._components: dict[int, DexpiObject] = {}
        self._templates: dict[int, DexpiObject] = {}
        self._streams: dict[int, DexpiObject] = {}
        self._identifiers: set[str] = {i.ID for i in self.units}

        self._note(model, 'system', {'ID': system.ID})
        for unit in self.units:
            self._export_unit(unit)

        root = DexpiObject(ENGINEERING_MODEL)
        root.add('ConceptualModel', model)
        timestamp = timestamp or datetime.now(timezone.utc)
        root.data.update({
            'ExportDateTime': timestamp.replace(microsecond=0),
            'OriginatingSystemName': 'BioSTEAM',
            'OriginatingSystemVendorName': 'BioSTEAM Development Group',
            'OriginatingSystemVersion': bst.__version__,
        })
        self.document = DexpiDocument(root, name=system.ID or 'BioSTEAM')

    # %% Annotations

    def _note(self, obj: DexpiObject, kind: str, payload: dict) -> None:
        note = create_note(kind, payload)
        self.model.add('Notes', note)
        obj.refer('ReferencedNotes', note)

    def _unique_identifier(self, base: str) -> str:
        identifier, count = base, 1
        while identifier in self._identifiers:
            count += 1
            identifier = f'{base}_{count}'
        self._identifiers.add(identifier)
        return identifier

    # %% Materials

    def _component(self, chemical: tmo.Chemical) -> DexpiObject:
        key = id(chemical)
        if key in self._components:
            return self._components[key]
        identity = ChemicalIdentity.from_chemical(chemical)
        kind = (
            'PureMaterialComponent' if identity.pure
            else 'CustomMaterialComponent'
        )
        component = DexpiObject(process_type(kind), data={
            'Identifier': identity.ID,
            'Label': identity.name or identity.ID,
        })
        if identity.pure:
            component.data['ChEBI_Identifier'] = identity.ChEBI
            component.data['IUPAC_Identifier'] = identity.IUPAC_name
        for name, context in _IDENTIFIER_CONTEXTS.items():
            value = getattr(identity, name)
            if value:
                component.add('PersistentIdentifiers', DexpiObject(
                    PERSISTENT_IDENTIFIER,
                    data={'Context': context, 'Value': str(value)},
                ))
        payload = {
            name: value for name, value in identity.to_dict().items()
            if name in ('formula', 'MW', 'phase', 'ChEBI')
        }
        self._note(component, 'chemical', payload)
        self.model.add('MaterialComponents', component)
        self._components[key] = component
        return component

    def _template(self, chemicals: tmo.Chemicals) -> DexpiObject:
        key = id(chemicals)
        if key in self._templates:
            return self._templates[key]
        components = DexpiObject(process_type('ListOfMaterialComponents'))
        components.refer('Component', *[self._component(i) for i in chemicals])
        self.model.add('ListsOfMaterialComponents', components)
        identifier = f'chemicals_{len(self._templates) + 1}'
        template = DexpiObject(process_type('MaterialTemplate'), data={
            'Description': multi_language_string(
                'Components of a BioSTEAM property package'
            ),
            'Identifier': identifier,
            'Label': identifier,
            'NumberOfMaterialComponents': len(chemicals),
            'NumberOfPhases': 1,
            'PhaseLabel': [],
        })
        template.refer('ListOfComponents', components)
        self.model.add('MaterialTemplates', template)
        self._templates[key] = template
        return template

    def _register_phases(self, template: DexpiObject, phases) -> None:
        labels = template.data['PhaseLabel']
        for phase in phases:
            if phase not in labels:
                labels.append(phase)
        template.data['NumberOfPhases'] = max(len(labels), 1)

    @staticmethod
    def _composition(mol) -> DexpiObject:
        composition = DexpiObject(process_type('Composition'), data={
            'Display': process_enumeration(
                'CompositionDisplay', 'AbsoluteValue'),
        })
        composition.add(
            'MoleFlow', *[_quantity(float(i), 'mole flow') for i in mol]
        )
        return composition

    def _state_type(self, identifier: str, label: str,
                    mol=None) -> DexpiObject:
        state = DexpiObject(process_type('MaterialStateType'), data={
            'Description': multi_language_string(f'Material state {label}'),
            'Identifier': identifier,
            'Label': label,
        })
        if mol is not None:
            state.add('Composition', self._composition(mol))
        return state

    def _material_state(self, stream) -> DexpiObject:
        ID = stream.ID
        state = DexpiObject(process_type('MaterialState'), data={
            'Description': f'Material state of stream {ID}',
            'Identifier': f'{ID}.state',
            'Label': ID,
        })
        phases = stream.phases
        if len(phases) == 1:
            state.add('State', self._state_type(
                f'{ID}.total', 'total', stream.mol))
            state.add('Phases', self._state_type(
                f'{ID}.{phases[0]}', phases[0]))
        else:
            imol = stream.imol
            state.add('State', self._state_type(
                f'{ID}.total', 'total', stream.mol))
            for phase in phases:
                state.add('Phases', self._state_type(
                    f'{ID}.{phase}', phase, imol[phase]
                ))
        self.model.add('MaterialStates', state)
        return state

    # %% Topology

    def _stream(self, stream) -> DexpiObject:
        key = id(stream)
        if key in self._streams:
            return self._streams[key]
        template = self._template(stream.chemicals)
        self._register_phases(template, stream.phases)
        obj = DexpiObject(process_type('Stream'), data={
            'Identifier': stream.ID,
            'Label': stream.ID,
        })
        obj.add('Temperature', _quantity(stream.T, 'temperature'))
        obj.add('Pressure', _quantity(stream.P, 'pressure'))
        obj.add('MassFlow', _quantity(stream.F_mass, 'mass flow'))
        obj.refer('MaterialStateReference', self._material_state(stream))
        obj.refer('MaterialTemplateReference', template)
        if stream.price:
            self._note(obj, 'stream', {'price': float(stream.price)})
        self.model.add('ProcessConnections', obj)
        self._streams[key] = obj
        self._identifiers.add(stream.ID)
        return obj

    def _port(self, step: DexpiObject, stream, direction,
              identifier: str) -> None:
        connection = self._stream(stream)
        port = DexpiObject(process_type('MaterialPort'), data={
            'Identifier': identifier,
            'NominalDirection': direction,
        })
        port.refer('ConnectorReference', connection)
        port.refer('MaterialTemplateReference',
                   self._template(stream.chemicals))
        step.add('Ports', port)
        connection.refer('Source' if direction is _OUTLET else 'Target', port)

    def _boundary(self, kind: str, stream) -> None:
        """Add the source or sink step that terminates a stream."""
        prefix = 'feed' if kind == 'Source' else 'product'
        identifier = self._unique_identifier(f'{prefix}.{stream.ID}')
        step = DexpiObject(process_type(kind), data={
            'Identifier': identifier,
            'Label': stream.ID,
        })
        if kind == 'Source':
            self._port(step, stream, _OUTLET, 'out0')
        else:
            self._port(step, stream, _INLET, 'in0')
        self.model.add('ProcessSteps', step)

    def _export_unit(self, unit) -> None:
        registration = self.registry.lookup_unit_class(type(unit))
        codec = registration.codec
        dexpi_class = codec.dexpi_class(unit, registration.dexpi_class)
        step = DexpiObject(process_type(dexpi_class), data={
            'Description': multi_language_string(unit.line),
            'Identifier': unit.ID,
            'Label': unit.ID,
        })
        for i, stream in enumerate(unit.ins):
            if not stream:
                continue  # Empty slot
            self._port(step, stream, _INLET, f'in{i}')
            if stream.source not in self._unit_set:
                self._boundary('Source', stream)
        for i, stream in enumerate(unit.outs):
            if not stream:
                continue
            self._port(step, stream, _OUTLET, f'out{i}')
            if stream.sink not in self._unit_set:
                self._boundary('Sink', stream)
        for name, (value, kind) in codec.step_quantities(unit).items():
            step.add(name, _quantity(value, kind))
        step.data.update(codec.step_data(unit))
        codec.annotate_step(unit, step)
        parameters, omitted = codec.encode_parameters(unit)
        payload = {'class': type(unit).__name__, 'parameters': parameters}
        if omitted:
            payload['omitted'] = omitted
        self._note(step, 'unit', payload)
        self.model.add('ProcessSteps', step)
