# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Import of DEXPI process models as BioSTEAM systems.

"""
from __future__ import annotations

import re
import warnings
from contextlib import contextmanager
from typing import Iterable, Optional, Union

import thermosteam as tmo

from ..chemicals import (
    ChemicalIdentity,
    ResolverChain,
    UserChemicalsResolver,
    default_resolver,
    is_cas_number,
    normalize_chebi,
)
from . import quantities
from .codecs import default_registry
from .generic import GenericProcessStep
from .model import DexpiDocument, DexpiObject, read_text
from .notes import read_note
from .registry import DexpiWarning, MappingError, Registration, UnitRegistry

__all__ = ('DexpiImportError', 'DexpiImporter')

_PORT_INDEX = re.compile(r'(\d+)$')
_BOUNDARY_CLASSES = frozenset({'Source', 'Sink'})
_KNOWN_PHASES = frozenset('slgLS')


class DexpiImportError(ValueError):
    """Raised when a document cannot be interpreted as a BioSTEAM system."""


@contextmanager
def _flowsheet(flowsheet):
    """Temporarily make `flowsheet` the main flowsheet."""
    import biosteam as bst
    if flowsheet is None:
        yield
        return
    main = bst.main_flowsheet
    previous = main.ID
    main.set_flowsheet(flowsheet)
    try:
        yield
    finally:
        main.set_flowsheet(previous)


def _active_thermo() -> Optional[tmo.Thermo]:
    try:
        return tmo.settings.get_thermo()
    except RuntimeError:
        return None


class DexpiImporter:
    """
    Translate a DEXPI document into a BioSTEAM system.

    Parameters
    ----------
    registry :
        Mapping between unit classes and DEXPI classes. Defaults to
        :data:`~biosteam.interop.dexpi.default_registry`.
    chemicals :
        Chemicals (or a thermodynamic property package) to use for the
        components of the document, matched by ID, CAS number or InChIKey.
        Components that are not covered are resolved with `resolver`.
    resolver :
        Resolver for the components of the document. Defaults to
        :func:`~biosteam.interop.chemicals.default_resolver`.
    flowsheet :
        Flowsheet (or the ID of a new one) in which units and streams are
        registered. Defaults to the main flowsheet.

    """

    def __init__(self, registry: Optional[UnitRegistry] = None,
                 chemicals: Union[
                     tmo.Thermo, Iterable[tmo.Chemical], None] = None,
                 resolver: Optional[ResolverChain] = None,
                 flowsheet=None):
        self.registry = registry or default_registry
        self.chemicals = chemicals
        self.resolver = resolver
        self.flowsheet = flowsheet

    def load(self, document: DexpiDocument, ID: Optional[str] = None):
        """
        Return the system described by `document`.

        Parameters
        ----------
        document :
            A document with a process model.
        ID :
            ID of the system. Defaults to the ID stored in the document.

        Raises
        ------
        DexpiImportError
            If the document has no process model or no process steps.

        """
        with _flowsheet(self.flowsheet):
            return _Import(self, document, ID).system


class _Import:
    """The state of a single import."""

    def __init__(self, importer: DexpiImporter, document: DexpiDocument, ID):
        import biosteam as bst
        self.registry = importer.registry
        model = document.process_model
        if model is None:
            raise DexpiImportError(
                "document does not contain a DEXPI process model"
            )
        self.model = model
        steps = [
            i for i in self._leaf_steps(model.children('ProcessSteps'))
            if i.class_name not in _BOUNDARY_CLASSES
        ]
        if not steps:
            raise DexpiImportError(
                "process model does not define any process steps")

        self._thermo: dict[int, tmo.Thermo] = {}
        self._default_thermo = self._load_thermo(importer)
        self._port_streams = self._index_ports()
        self._streams: dict[int, tmo.Stream] = {}

        units = [self._create_unit(i) for i in steps]
        if ID is None:
            ID = (read_note(model) or {}).get('ID') or document.name or None
        self.system = bst.System.from_units(ID, units)

    # %% Thermodynamic property package

    @staticmethod
    def _identity(component: DexpiObject) -> ChemicalIdentity:
        note = read_note(component) or {}
        identifiers = {
            i.get('Context'): i.get('Value')
            for i in component.children('PersistentIdentifiers')
        }
        ID = component.get('Identifier') or read_text(component.get('Label'))
        if not ID:
            raise DexpiImportError("material component without an identifier")
        label = component.get('Label')
        CAS = identifiers.get('CAS')
        if CAS is None and is_cas_number(ID):
            CAS = ID
        return ChemicalIdentity(
            ID=ID,
            name=label if isinstance(label, str) and label != ID else None,
            CAS=CAS,
            InChIKey=identifiers.get('InChIKey'),
            InChI=identifiers.get('InChI'),
            SMILES=identifiers.get('SMILES'),
            PubChem=identifiers.get('PubChem'),
            ChEBI=normalize_chebi(
                component.get('ChEBI_Identifier') or note.get('ChEBI')
            ),
            IUPAC_name=component.get('IUPAC_Identifier'),
            formula=note.get('formula'),
            MW=note.get('MW'),
            phase=note.get('phase'),
            pure=component.class_name != 'CustomMaterialComponent',
        )

    def _template_components(self, template: DexpiObject) -> list[DexpiObject]:
        components = template.target('ListOfComponents')
        return components.targets('Component') if components else []

    def _load_thermo(self, importer: DexpiImporter) -> tmo.Thermo:
        """Create the property packages of all material templates."""
        templates = self.model.children('MaterialTemplates')
        if not templates:
            # Without templates, all material components form one package.
            components = self.model.children('MaterialComponents')
            return self._create_thermo(importer, components)
        for template in templates:
            self._thermo[id(template)] = self._create_thermo(
                importer, self._template_components(template)
            )
        return self._thermo[id(templates[0])]

    def _create_thermo(self, importer: DexpiImporter,
                       components) -> tmo.Thermo:
        identities = [self._identity(i) for i in components]
        IDs = [i.ID for i in identities]
        if len(set(IDs)) != len(IDs):
            raise DexpiImportError(
                "material components must have unique identifiers")
        given = importer.chemicals
        if isinstance(given, tmo.Thermo) and given.chemicals.IDs == tuple(IDs):
            return given
        active = _active_thermo()
        if given is None and importer.resolver is None:
            if active is not None and active.chemicals.IDs == tuple(IDs):
                return active
            if not IDs:
                if active is None:
                    raise DexpiImportError(
                        "document does not define material components; set "
                        "a property package with settings.set_thermo first"
                    )
                return active
        resolver = importer.resolver
        if resolver is None:
            if isinstance(given, tmo.Thermo):
                given = given.chemicals
            resolver = default_resolver(chemicals=given)
        elif given is not None:
            if isinstance(given, tmo.Thermo):
                given = given.chemicals
            resolver = ResolverChain([UserChemicalsResolver(given), resolver])
        thermo = tmo.Thermo(resolver.create_chemicals(identities))
        if active is None:
            tmo.settings.set_thermo(thermo)
        return thermo

    def _thermo_of(self, obj: DexpiObject) -> tmo.Thermo:
        template = obj.target('MaterialTemplateReference')
        if template is None:
            return self._default_thermo
        return self._thermo.get(id(template), self._default_thermo)

    # %% Streams

    def _index_ports(self) -> dict[int, DexpiObject]:
        """Map each port to the stream connected to it."""
        port_streams = {}
        for connection in self.model.children('ProcessConnections'):
            if connection.class_name != 'Stream':
                continue
            for name in ('Source', 'Target'):
                port = connection.target(name)
                if port is not None:
                    port_streams[id(port)] = connection
        return port_streams

    def _connection(self, port: DexpiObject) -> Optional[DexpiObject]:
        connection = self._port_streams.get(id(port))
        if connection is None:
            connection = port.target('ConnectorReference')
            if connection is not None and connection.class_name != 'Stream':
                return None  # Energy and information flows are not imported.
        return connection

    def _flows(self, state: DexpiObject, size: int) -> Optional[list[float]]:
        """Return the component mole flows of a material state [kmol/hr]."""
        composition = state.child('Composition') if state else None
        if composition is None:
            return None
        flows = [
            quantities.from_dexpi(i, 'mole flow') or 0.
            for i in composition.children('MoleFlow')
        ]
        if not flows:
            return None
        if len(flows) != size:
            raise DexpiImportError(
                f"composition of '{state.get('Identifier')}' has {len(flows)} "
                f"values, but its material template has {size} components"
            )
        return flows

    def _stream(self, connection: DexpiObject):
        import biosteam as bst
        key = id(connection)
        if key in self._streams:
            return self._streams[key]
        thermo = self._thermo_of(connection)
        size = thermo.chemicals.size
        ID = (
            connection.get('Identifier')
            or read_text(connection.get('Label')) or ''
        )
        conditions = {}
        T = quantities.first(connection, 'Temperature', 'temperature')
        if T is not None:
            conditions['T'] = T
        P = quantities.first(connection, 'Pressure', 'pressure')
        if P is not None:
            conditions['P'] = P

        state = connection.target('MaterialStateReference')
        total = state.child('State') if state else None
        phases = [
            i for i in (state.children('Phases') if state else [])
            if i.get('Label') in _KNOWN_PHASES
        ]
        phase_flows = [(i.get('Label'), self._flows(i, size)) for i in phases]
        if len(phases) > 1 and all(flows for _, flows in phase_flows):
            stream = bst.MultiStream(
                ID, phases=[i for i, _ in phase_flows], thermo=thermo,
                **conditions
            )
            for phase, flows in phase_flows:
                stream.imol[phase] = flows
        else:
            flows = self._flows(total, size)
            phase = 'l'
            if phases:
                phase = phases[0].get('Label')
                if flows is None:
                    flows = phase_flows[0][1]
            stream = bst.Stream(ID, phase=phase, thermo=thermo, **conditions)
            if flows is not None:
                stream.mol[:] = flows
            else:
                self._set_total_flow(stream, connection)
        note = read_note(connection)
        if note and note.get('price'):
            stream.price = float(note['price'])
        self._streams[key] = stream
        return stream

    @staticmethod
    def _set_total_flow(stream, connection: DexpiObject) -> None:
        """Warn about a stream whose composition is not defined."""
        mass_flow = quantities.first(connection, 'MassFlow', 'mass flow')
        if mass_flow:
            warnings.warn(
                f"stream '{stream.ID}' defines a mass flow but no "
                "composition; it is imported empty",
                DexpiWarning, stacklevel=6,
            )

    # %% Unit operations

    def _leaf_steps(self, steps: Iterable[DexpiObject]):
        """Iterate over process steps, replacing steps by their sub-steps."""
        for step in steps:
            children = step.children('SubProcessSteps')
            if children:
                yield from self._leaf_steps(children)
            else:
                yield step

    def _connected_streams(self, step: DexpiObject) -> tuple[list, list]:
        """Return the inlet and outlet streams of a step by position."""
        ins, outs = {}, {}
        for port in step.children('Ports'):
            connection = self._connection(port)
            if connection is None:
                continue
            direction = port.get('NominalDirection')
            is_outlet = getattr(direction, 'literal', None) == 'Outlet'
            if direction is None:
                is_outlet = connection.target('Source') is port
            slots = outs if is_outlet else ins
            match = _PORT_INDEX.search(port.get('Identifier') or '')
            index = int(match.group(1)) if match else len(slots)
            while index in slots:
                index += 1
            slots[index] = self._stream(connection)

        def ordered(slots):
            if not slots:
                return []
            if max(slots) >= 2 * len(slots) + 8:
                # Identifiers are not positions; keep the order only.
                return [slots[i] for i in sorted(slots)]
            return [slots.get(i) for i in range(max(slots) + 1)]

        return ordered(ins), ordered(outs)

    def _registration(self, step: DexpiObject, note: Optional[dict]
                      ) -> tuple[Optional[Registration], bool]:
        """Return the registration for a step and whether its annotation
        applies to it."""
        if note and note.get('class'):
            unit_class = self.registry.resolve_class_name(str(note['class']))
            if unit_class is not None:
                codec = self.registry.lookup_unit_class(unit_class).codec
                return Registration(unit_class, step.class_name, codec), True
            warnings.warn(
                f"unit class '{note['class']}' of step "
                f"'{step.get('Identifier')}' is not available; using the "
                f"default model for DEXPI class '{step.class_name}'",
                DexpiWarning, stacklevel=5,
            )
        return self.registry.lookup_dexpi_class(step.class_name), False

    def _create_unit(self, step: DexpiObject):
        ID = step.get('Identifier') or read_text(step.get('Label')) or ''
        ins, outs = self._connected_streams(step)
        streams = [i for i in (*ins, *outs) if i is not None]
        thermo = streams[0].thermo if streams else self._default_thermo
        note = read_note(step)
        if note and note.get('kind') not in (None, 'unit'):
            note = None
        registration, annotated = self._registration(step, note)
        if registration is not None:
            codec = registration.codec
            try:
                if annotated:
                    kwargs = codec.decode_parameters(
                        note.get('parameters') or {}, thermo.chemicals
                    )
                else:
                    kwargs = codec.parameters_from_step(step, thermo.chemicals)
                unit = registration.unit_class(
                    ID, ins=ins, outs=outs, thermo=thermo, **kwargs
                )
            except (MappingError, TypeError, ValueError, KeyError) as error:
                warnings.warn(
                    f"step '{ID}' could not be created as a "
                    f"'{registration.unit_class.__name__}' unit ({error}); "
                    "using a generic process step instead",
                    DexpiWarning, stacklevel=4,
                )
            else:
                codec.finalize(unit, step)
                if isinstance(unit, GenericProcessStep):
                    unit.dexpi_class = step.class_name
                    unit.dexpi_step = step
                return unit
        unit = GenericProcessStep(
            ID, ins=ins, outs=outs, thermo=thermo, dexpi_class=step.class_name,
        )
        unit.dexpi_step = step
        return unit
