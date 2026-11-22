# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Codecs for BioSTEAM's built-in unit operations and the default registry.

"""
from __future__ import annotations

import copy
import warnings

from ... import units
from ..reactions import iter_rates
from . import quantities
from .generic import GenericProcessStep
from .model import DexpiObject, process_enumeration
from .registry import DexpiWarning, MappingError, UnitCodec, UnitRegistry

__all__ = (
    'FlashCodec',
    'GenericStepCodec',
    'HeatExchangerCodec',
    'PumpCodec',
    'ReactorCodec',
    'SplitterCodec',
    'ValveCodec',
    'create_default_registry',
    'default_registry',
)

#: Attributes of a process step that the exporter always writes itself.
_STRUCTURAL_DATA = frozenset({'Identifier', 'Label'})
_STRUCTURAL_COMPONENTS = frozenset({'Ports', 'SubProcessSteps'})


def _operating_conditions(step: DexpiObject) -> dict:
    """Return the temperature and pressure of a step as unit parameters."""
    conditions = {}
    T = quantities.first(step, 'Temperature', 'temperature')
    if T is not None:
        conditions['T'] = T
    P = quantities.first(step, 'Pressure', 'pressure')
    if P is not None:
        conditions['P'] = P
    return conditions


class ReactorCodec(UnitCodec):
    """Codec for stirred tank reactors (``ReactingChemicals`` steps)."""

    #: Residence time assumed when a document does not define one [hr].
    default_residence_time = 1.

    #: Flow pattern assumed for each DEXPI reaction process type.
    flow_patterns = {
        'Tank': 'mixed',
        'Tubular': 'plug',
        'PackedBed': 'plug',
    }

    def step_data(self, unit):
        plug_flow = getattr(unit, 'flow_pattern', None) == 'plug'
        method = 'Tubular' if plug_flow else 'Tank'
        return {'Method': process_enumeration('ReactionProcessType', method)}

    def step_quantities(self, unit):
        step_quantities = super().step_quantities(unit)
        volume = (unit.design_results.get('Reactor volume')
                  or getattr(unit, 'dexpi_volume', None))
        if volume:
            step_quantities['Volume'] = (float(volume), 'volume')
        return step_quantities

    def parameters_from_step(self, step, chemicals):
        warnings.warn(
            f"reactor '{step.get('Identifier')}' has no BioSTEAM annotation; "
            "it is created without reactions and with a residence time of "
            f"{self.default_residence_time:g} hr",
            DexpiWarning, stacklevel=5,
        )
        return {
            **_operating_conditions(step),
            'tau': self.default_residence_time,
        }

    def finalize(self, unit, step=None):
        # Kinetic rates depend on the reactor volume; let them read the
        # residence time of the reactor they belong to.
        for rate in iter_rates(getattr(unit, 'reactions', None)):
            bind = getattr(rate, 'bind', None)
            if callable(bind):
                bind(unit)
        if step is None:
            return
        # Keep what the document says about the reactor but BioSTEAM has
        # no parameter for: the flow pattern, which selects the kinetic
        # model, and the reactor volume.
        method = getattr(step.get('Method'), 'literal', None)
        unit.flow_pattern = self.flow_patterns.get(method, 'mixed')
        unit.dexpi_volume = quantities.first(step, 'Volume', 'volume')


class SplitterCodec(UnitCodec):
    """Codec for splitters (``SplittingMaterial`` steps)."""

    excluded_parameters = frozenset({'order'})

    def read_parameter(self, unit, name):
        if name == 'split':
            return unit.isplit
        return super().read_parameter(unit, name)

    def parameters_from_step(self, step, chemicals):
        warnings.warn(
            f"splitter '{step.get('Identifier')}' has no BioSTEAM annotation; "
            "assuming an even split",
            DexpiWarning, stacklevel=5,
        )
        return {'split': 0.5}


class PumpCodec(UnitCodec):
    """Codec for pumps (``Pumping`` steps)."""

    def parameters_from_step(self, step, chemicals):
        conditions = _operating_conditions(step)
        return {'P': conditions['P']} if 'P' in conditions else {}


class HeatExchangerCodec(UnitCodec):
    """Codec for utility heat exchangers (``ExchangingThermalEnergy``)."""

    def step_data(self, unit):
        return {'Method': process_enumeration('HeatExchangeMethod', 'Generic')}

    def parameters_from_step(self, step, chemicals):
        conditions = _operating_conditions(step)
        if 'T' not in conditions:
            raise MappingError(
                f"heat exchanger '{step.get('Identifier')}' does not define "
                "an outlet temperature"
            )
        return {'T': conditions['T']}


class FlashCodec(UnitCodec):
    """Codec for flash vessels (``SeparatingByFlash`` steps)."""

    def parameters_from_step(self, step, chemicals):
        conditions = _operating_conditions(step)
        conditions.setdefault('P', 101325.)
        if 'T' not in conditions:
            conditions['Q'] = 0.  # Adiabatic unless a temperature is given.
        return conditions


class ValveCodec(UnitCodec):
    """Codec for isenthalpic valves (``RegulatingFlow`` steps)."""

    def parameters_from_step(self, step, chemicals):
        conditions = _operating_conditions(step)
        if 'P' not in conditions:
            raise MappingError(
                f"valve '{step.get('Identifier')}' does not define an "
                "outlet pressure"
            )
        return {'P': conditions['P']}


class GenericStepCodec(UnitCodec):
    """
    Codec for :class:`GenericProcessStep`, which restores the DEXPI class
    and attributes of the process step the unit was imported from.
    """

    def dexpi_class(self, unit, default):
        return unit.dexpi_class

    def step_quantities(self, unit):
        return {}

    def annotate_step(self, unit, step):
        original = unit.dexpi_step
        if original is None:
            return
        for name, value in original.data.items():
            if name not in _STRUCTURAL_DATA:
                step.data.setdefault(name, copy.deepcopy(value))
        for name, children in original.components.items():
            if name in _STRUCTURAL_COMPONENTS or name in step.components:
                continue
            step.components[name] = [_copy_without_ids(i) for i in children]


def _copy_without_ids(obj: DexpiObject) -> DexpiObject:
    """Return a copy of an object tree with fresh identifiers.

    References are dropped, since their targets belong to the document the
    object was read from.
    """
    return DexpiObject(
        obj.type,
        data=copy.deepcopy(obj.data),
        components={
            name: [_copy_without_ids(i) for i in children]
            for name, children in obj.components.items()
        },
    )


def create_default_registry() -> UnitRegistry:
    """Return a new registry with the mappings of the built-in units."""
    registry = UnitRegistry()
    register = registry.register
    register(units.Mixer, 'MixingSimple', import_default_for=('Mixing',))
    register(units.Splitter, 'SplittingMaterial', SplitterCodec(),
             import_default=True)
    register(units.Pump, 'Pumping', PumpCodec(), import_default=True)
    register(units.HXutility, 'ExchangingThermalEnergy', HeatExchangerCodec(),
             import_default=True,
             import_default_for=(
                 'SupplyingThermalEnergy', 'RemovingThermalEnergy'))
    register(units.Flash, 'SeparatingByFlash', FlashCodec(),
             import_default=True)
    register(units.StorageTank, 'StoringInTank', import_default=True)
    register(units.MixTank, 'StoringInTank')
    register(units.IsenthalpicValve, 'RegulatingFlow', ValveCodec(),
             import_default=True)
    reactor_codec = ReactorCodec()
    register(units.AbstractStirredTankReactor, 'ReactingChemicals',
             reactor_codec)
    register(units.SinglePhaseReactor, 'ReactingChemicals', reactor_codec,
             import_default=True)
    register(GenericProcessStep, 'ProcessStep', GenericStepCodec())
    return registry


#: Registry used by default for export and import.
default_registry = create_default_registry()
