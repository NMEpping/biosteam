# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Conversion between DEXPI physical quantities and BioSTEAM base units.

Each :class:`QuantityKind` ties one DEXPI unit enumeration to the unit
BioSTEAM uses internally for that kind of quantity, e.g. temperatures are
exchanged with ``TemperatureUnit`` literals and held in Kelvin.

"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from .model import Aggregate, DexpiObject, physical_quantity, read_quantity

__all__ = (
    'QuantityKind', 'UnitConversionError', 'KINDS', 'from_dexpi', 'to_dexpi',
)


class UnitConversionError(ValueError):
    """Raised when a DEXPI unit is unknown or does not fit a quantity."""


@dataclass(frozen=True)
class QuantityKind:
    """
    A kind of physical quantity exchanged with DEXPI.

    Parameters
    ----------
    enumeration :
        Name of the DEXPI unit enumeration.
    base_unit :
        Unit used by BioSTEAM for this kind of quantity.
    literal :
        DEXPI unit literal written on export.
    factors :
        For each DEXPI unit literal, the scale and offset such that
        ``base = value * scale + offset``.

    """

    enumeration: str
    base_unit: str
    literal: str
    factors: dict

    def to_base(self, value: float, literal: Optional[str]) -> float:
        """Convert a value in DEXPI units to BioSTEAM base units."""
        if literal is None:
            return value
        try:
            scale, offset = self.factors[literal]
        except KeyError:
            raise UnitConversionError(
                f"unit '{literal}' is not a known {self.enumeration}"
            ) from None
        return value * scale + offset

    def from_base(self, value: float) -> float:
        """Convert a value in BioSTEAM base units to the export unit."""
        scale, offset = self.factors[self.literal]
        return (value - offset) / scale


def _kind(enumeration, base_unit, literal, **factors):
    factors = {
        name: factor if isinstance(factor, tuple) else (factor, 0.)
        for name, factor in factors.items()
    }
    return QuantityKind(enumeration, base_unit, literal, factors)


_LB = 0.45359237  # kg
_FT3 = 0.028316846592  # m3

#: Quantity kinds by name.
KINDS = {
    'temperature': _kind(
        'TemperatureUnit', 'K', 'Kelvin',
        Kelvin=1.,
        DegreeCelsius=(1., 273.15),
        DegreeFahrenheit=(5. / 9., 459.67 * 5. / 9.),
    ),
    'pressure': _kind(
        'PressureAbsoluteUnit', 'Pa', 'Pascal',
        Pascal=1.,
        Kilopascal=1e3,
        Megapascal=1e6,
        Bar=1e5,
        Millibar=1e2,
        PoundForcePerInchSquared=6894.757293168361,
    ),
    'mass flow': _kind(
        'MassFlowRateUnit', 'kg/hr', 'KilogramPerHour',
        KilogramPerHour=1.,
        KilogramPerMinute=60.,
        KilogramPerSecond=3600.,
        PoundMassPerHour=_LB,
        PoundMassPerMinute=_LB * 60.,
        PoundMassPerSecond=_LB * 3600.,
    ),
    'mole flow': _kind(
        'MoleFlowRateUnit', 'kmol/hr', 'KilomolePerSecond',
        KilomolePerSecond=3600.,
        PoundMolePerSecond=_LB * 3600.,
    ),
    'volume flow': _kind(
        'VolumeFlowRateUnit', 'm3/hr', 'MetreCubedPerHour',
        MetreCubedPerHour=1.,
        MetreCubedPerMinute=60.,
        MetreCubedPerSecond=3600.,
        MetreCubedPerDay=1. / 24.,
        LitrePerSecond=3.6,
        FootCubedPerHour=_FT3,
        FootCubedPerMinute=_FT3 * 60.,
    ),
    'volume': _kind(
        'VolumeUnit', 'm3', 'MetreCubed',
        MetreCubed=1.,
        DecimetreCubed=1e-3,
        Litre=1e-3,
        CentimetreCubed=1e-6,
        FootCubed=_FT3,
        UsGallon=3.785411784e-3,
        UsFluidOunce=2.95735295625e-5,
    ),
    'power': _kind(
        'PowerUnit', 'kW', 'Kilowatt',
        Kilowatt=1.,
        Watt=1e-3,
        Megawatt=1e3,
    ),
    'ratio': _kind(
        'RatioUnit', '', 'Fraction',
        Fraction=1.,
        Percent=1e-2,
        PartsPerThousand=1e-3,
        PartsPerMillion=1e-6,
        PartsPerBillion=1e-9,
    ),
}


def to_dexpi(value: float, kind: str) -> Aggregate:
    """
    Return a DEXPI physical quantity for a value in BioSTEAM base units.

    Examples
    --------
    >>> from biosteam.interop.dexpi import quantities
    >>> quantity = quantities.to_dexpi(7200., 'mole flow')
    >>> quantity.data['Value'], quantity.data['Unit'].literal
    (2.0, 'KilomolePerSecond')

    """
    quantity_kind = KINDS[kind]
    return physical_quantity(
        quantity_kind.from_base(float(value)),
        quantity_kind.enumeration,
        quantity_kind.literal,
    )


def from_dexpi(value: Any, kind: str) -> Optional[float]:
    """
    Return a quantity in BioSTEAM base units, or None if `value` holds
    no number.

    Parameters
    ----------
    value :
        A qualified value object, a physical quantity or a plain number.
    kind :
        Name of the quantity kind (a key of :data:`KINDS`).

    Examples
    --------
    >>> from biosteam.interop.dexpi import model, quantities
    >>> celsius = model.physical_quantity(
    ...     25., 'TemperatureUnit', 'DegreeCelsius')
    >>> quantities.from_dexpi(celsius, 'temperature')
    298.15

    """
    quantity = read_quantity(value)
    if quantity is None:
        return None
    number, unit = quantity
    return KINDS[kind].to_base(number, unit.literal if unit else None)


def first(obj: DexpiObject, name: str, kind: str) -> Optional[float]:
    """Return the first value of a quantity attribute of an object."""
    for child in obj.children(name):
        value = from_dexpi(child, kind)
        if value is not None:
            return value
    return None
