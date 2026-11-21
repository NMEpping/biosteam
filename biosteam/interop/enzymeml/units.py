# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Conversion of EnzymeML units to the unit system of BioSTEAM.

EnzymeML describes a unit as a product of base units, each with a kind, an
exponent, a multiplier and a decimal scale (as in SBML). Kinetic parameters
are converted to one coherent system of units (kmol, m, hr, kg, K), so that
a rate law evaluated with concentrations in kmol/m3 returns a rate in
kmol/m3/hr.

"""
from __future__ import annotations

from typing import Any, Mapping, Optional

__all__ = ('UnitError', 'conversion_factor', 'to_kelvin', 'unit_name')

#: Factor from each base unit kind to the coherent system (kmol, m, hr, kg, K).
_FACTORS = {
    'mole': 1e-3,
    'litre': 1e-3,
    'liter': 1e-3,
    'metre': 1.,
    'meter': 1.,
    'second': 1. / 3600.,
    'gram': 1e-3,
    'kilogram': 1.,
    'kelvin': 1.,
    'dimensionless': 1.,
    'item': 1.,
    'katal': 3.6,
}

#: Multipliers of the minute, hour and day in seconds.
_TIME_MULTIPLIERS = (60., 3600., 86400.)


class UnitError(ValueError):
    """Raised when an EnzymeML unit cannot be converted."""


def _kind(base_unit: Mapping[str, Any]) -> str:
    kind = base_unit.get('kind')
    kind = getattr(kind, 'value', kind)
    return str(kind).lower()


def unit_name(unit: Optional[Mapping[str, Any]]) -> Optional[str]:
    """Return the display name of an EnzymeML unit definition."""
    if not unit:
        return None
    return unit.get('name') or unit.get('id')


def conversion_factor(unit: Optional[Mapping[str, Any]]) -> float:
    """
    Return the factor that converts a value in an EnzymeML unit to the
    coherent system of units (kmol, m, hr, kg, K).

    Parameters
    ----------
    unit :
        An EnzymeML unit definition with a `base_units` list. A missing
        unit is treated as dimensionless.

    Examples
    --------
    One millimole per litre is 1e-3 kmol/m3:

    >>> from biosteam.interop.enzymeml.units import conversion_factor
    >>> mM = {'base_units': [
    ...     {'kind': 'mole', 'exponent': 1, 'multiplier': 1.0, 'scale': -3},
    ...     {'kind': 'litre', 'exponent': -1, 'multiplier': 1.0, 'scale': 0},
    ... ]}
    >>> round(conversion_factor(mM), 12)
    0.001

    One per second is 3600 per hour:

    >>> per_second = {'base_units': [
    ...     {'kind': 'second', 'exponent': -1, 'multiplier': 1.0, 'scale': 0},
    ... ]}
    >>> round(conversion_factor(per_second), 9)
    3600.0

    """
    if not unit:
        return 1.
    factor = 1.
    for base_unit in unit.get('base_units') or ():
        kind = _kind(base_unit)
        try:
            base = _FACTORS[kind]
        except KeyError:
            raise UnitError(
                f"base unit '{kind}' of unit '{unit_name(unit)}' is not "
                "supported"
            ) from None
        exponent = float(base_unit.get('exponent', 1) or 1)
        multiplier = float(base_unit.get('multiplier', 1.) or 1.)
        scale = float(base_unit.get('scale', 0.) or 0.)
        if (kind == 'second' and scale == 1.
                and multiplier in _TIME_MULTIPLIERS):
            # Some EnzymeML writers encode minutes and hours with a scale
            # of 1 in addition to the multiplier; the multiplier alone is
            # the number of seconds.
            scale = 0.
        factor *= (multiplier * 10. ** scale * base) ** exponent
    return factor


def to_kelvin(value: Optional[float],
              unit: Optional[Mapping[str, Any]]) -> Optional[float]:
    """
    Return a temperature in Kelvin.

    Examples
    --------
    >>> from biosteam.interop.enzymeml.units import to_kelvin
    >>> to_kelvin(30., {'base_units': [{'kind': 'celsius', 'exponent': 1}]})
    303.15

    """
    if value is None:
        return None
    kinds = [_kind(i) for i in (unit or {}).get('base_units') or ()]
    if kinds == ['celsius']:
        return value + 273.15
    return float(value)
