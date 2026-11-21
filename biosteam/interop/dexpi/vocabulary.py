# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
The class hierarchy of process steps in the DEXPI process package.

The hierarchy lets BioSTEAM interpret process steps it has no specific
mapping for through their more general parent class, e.g. a ``RotaryMixing``
step as ``Mixing``.

The class names and their relations are taken from the DEXPI Specification
2.0.1, (c) DEXPI e.V., licensed under the Creative Commons Attribution 4.0
International License (https://creativecommons.org/licenses/by/4.0/).

"""
from __future__ import annotations

from typing import Iterator

__all__ = ('PROCESS_STEP_PARENTS', 'ancestors', 'is_process_step')

#: Name of the root class of all process steps.
PROCESS_STEP = 'ProcessStep'

#: Parent class of each process step class.
PROCESS_STEP_PARENTS = {
    'Absorbing': 'SeparatingByPhysicalProcess',
    'Adsorbing': 'SeparatingByPhysicalProcess',
    'Agglomerating': 'IncreasingParticleSize',
    'BlowingDown': 'SteeringFlow',
    'Boiling': 'SupplyingThermalEnergy',
    'Coalescing': 'IncreasingParticleSize',
    'Compressing': 'GeneratingFlow',
    'Cooling': 'RemovingThermalEnergy',
    'Crushing': 'ReducingParticleSize',
    'Crystallizing': 'IncreasingParticleSize',
    'Cutting': 'ReducingParticleSize',
    'Distilling': 'SeparatingByThermalProcess',
    'Draining': 'SteeringFlow',
    'DrivingByEngine': 'SupplyingMechanicalEnergy',
    'DrivingByMotor': 'SupplyingMechanicalEnergy',
    'DrivingByTurbine': 'SupplyingMechanicalEnergy',
    'Drying': 'SeparatingByThermalProcess',
    'Emitting': 'ProcessStep',
    'Evaporating': 'SeparatingByThermalProcess',
    'ExchangingThermalEnergy': 'ProcessStep',
    'Extruding': 'FormingSolidMaterial',
    'FeedingMaterial': 'SteeringFlow',
    'Filtering': 'SeparatingMechanically',
    'Flaring': 'ProcessStep',
    'Flocculating': 'IncreasingParticleSize',
    'FormingSolidMaterial': 'ProcessStep',
    'GeneratingACPower': 'SupplyingElectricalEnergy',
    'GeneratingDCPower': 'SupplyingElectricalEnergy',
    'GeneratingFlow': 'ProcessStep',
    'GeneratingInFuelCell': 'SupplyingElectricalEnergy',
    'GeneratingSteam': 'SupplyingThermalEnergy',
    'Grinding': 'ReducingParticleSize',
    'HeatingElectrical': 'SupplyingThermalEnergy',
    'HeatingInFurnace': 'SupplyingThermalEnergy',
    'Humidifying': 'Mixing',
    'IncreasingParticleSize': 'ProcessStep',
    'Kneading': 'Mixing',
    'LimitingFlow': 'SteeringFlow',
    'Milling': 'ReducingParticleSize',
    'Mixing': 'ProcessStep',
    'MixingSimple': 'Mixing',
    'Packaging': 'ProcessStep',
    'Pelletizing': 'FormingSolidMaterial',
    'PreventingBackflow': 'SteeringFlow',
    'Pumping': 'GeneratingFlow',
    'ReactingChemicals': 'ProcessStep',
    'ReducingParticleSize': 'ProcessStep',
    'RegulatingFlow': 'SteeringFlow',
    'RelievingOverpressure': 'SteeringFlow',
    'RelievingVacuum': 'SteeringFlow',
    'RelievingVacuumAndOverpressure': 'SteeringFlow',
    'RemovingThermalEnergy': 'ProcessStep',
    'RotaryMixing': 'Mixing',
    'Separating': 'ProcessStep',
    'SeparatingByCentrifugalForce': 'SeparatingByPhaseSeparation',
    'SeparatingByContact': 'SeparatingByPhysicalProcess',
    'SeparatingByCyclonicMotion': 'SeparatingByPhaseSeparation',
    'SeparatingByElectromagneticForce': 'Separating',
    'SeparatingByElectrostaticForce': 'SeparatingByElectromagneticForce',
    'SeparatingByFlash': 'SeparatingByThermalProcess',
    'SeparatingByGravity': 'SeparatingByPhaseSeparation',
    'SeparatingByIonExchange': 'SeparatingByPhysicalProcess',
    'SeparatingByMagneticForce': 'SeparatingByElectromagneticForce',
    'SeparatingByPhaseSeparation': 'Separating',
    'SeparatingByPhysicalProcess': 'Separating',
    'SeparatingBySurfaceTension': 'SeparatingByPhysicalProcess',
    'SeparatingByThermalProcess': 'Separating',
    'SeparatingMechanically': 'Separating',
    'ShuttingOffFlow': 'SteeringFlow',
    'Sieving': 'SeparatingMechanically',
    'Sink': 'ProcessStep',
    'Skimming': 'SeparatingMechanically',
    'Source': 'ProcessStep',
    'Splitting': 'ProcessStep',
    'SplittingEnergy': 'Splitting',
    'SplittingMaterial': 'Splitting',
    'StabilizingDistilling': 'Distilling',
    'StaticMixing': 'Mixing',
    'SteeringFlow': 'ProcessStep',
    'StoringElectricalEnergy': 'StoringEnergy',
    'StoringEnergy': 'ProcessStep',
    'StoringFluids': 'StoringMaterial',
    'StoringInBattery': 'StoringElectricalEnergy',
    'StoringInPressureVessel': 'StoringFluids',
    'StoringInSilo': 'StoringSolids',
    'StoringInTank': 'StoringFluids',
    'StoringMaterial': 'ProcessStep',
    'StoringMechanicalEnergy': 'StoringEnergy',
    'StoringSolids': 'StoringMaterial',
    'StoringThermalEnergy': 'StoringEnergy',
    'StrippingDistilling': 'Distilling',
    'SupplyingElectricalEnergy': 'SupplyingEnergy',
    'SupplyingEnergy': 'ProcessStep',
    'SupplyingFluids': 'SupplyingMaterial',
    'SupplyingMaterial': 'ProcessStep',
    'SupplyingMechanicalEnergy': 'SupplyingEnergy',
    'SupplyingSolids': 'SupplyingMaterial',
    'SupplyingThermalEnergy': 'ProcessStep',
    'TransportingElectricalEnergy': 'TransportingEnergy',
    'TransportingEnergy': 'ProcessStep',
    'TransportingFluids': 'TransportingMaterial',
    'TransportingFluidsInChannel': 'TransportingFluids',
    'TransportingFluidsInHose': 'TransportingFluids',
    'TransportingFluidsInPipe': 'TransportingFluids',
    'TransportingMaterial': 'ProcessStep',
    'TransportingMechanicalEnergy': 'TransportingEnergy',
    'TransportingSolids': 'TransportingMaterial',
    'TransportingSolidsContinuously': 'TransportingSolids',
    'TransportingSolidsDiscontinuously': 'TransportingSolids',
    'TransportingThermalEnergy': 'TransportingEnergy',
    'VacuumDistilling': 'Distilling',
}


def is_process_step(class_name: str) -> bool:
    """
    Return whether a DEXPI class is a process step.

    Examples
    --------
    >>> from biosteam.interop.dexpi.vocabulary import is_process_step
    >>> is_process_step('Pumping'), is_process_step('Stream')
    (True, False)

    """
    return class_name == PROCESS_STEP or class_name in PROCESS_STEP_PARENTS


def ancestors(class_name: str) -> Iterator[str]:
    """
    Iterate over a process step class and its parents, most specific first.

    Classes that are not part of the DEXPI vocabulary are treated as direct
    children of ``ProcessStep``.

    Examples
    --------
    >>> from biosteam.interop.dexpi.vocabulary import ancestors
    >>> list(ancestors('RotaryMixing'))
    ['RotaryMixing', 'Mixing', 'ProcessStep']

    """
    seen = set()
    while class_name and class_name not in seen:
        yield class_name
        seen.add(class_name)
        if class_name == PROCESS_STEP:
            return
        class_name = PROCESS_STEP_PARENTS.get(class_name, PROCESS_STEP)
