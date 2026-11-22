# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
The mapping between BioSTEAM unit operations and DEXPI process steps.

A :class:`UnitRegistry` records which DEXPI class represents each BioSTEAM
unit class and which BioSTEAM class models each DEXPI class. How a unit is
translated is delegated to a :class:`UnitCodec`, so that new unit
operations can be supported by registering a codec, without changing the
exporter or importer.

"""
from __future__ import annotations

import inspect
import warnings
from dataclasses import dataclass
from numbers import Real
from typing import Any, Mapping, Optional

import thermosteam as tmo

from ..._unit import Unit
from .model import DexpiObject
from .values import ValueNotSerializable, decode_value, encode_value
from .vocabulary import PROCESS_STEP, ancestors

__all__ = (
    'DexpiWarning',
    'MappingError',
    'Registration',
    'UnitCodec',
    'UnitRegistry',
)


class DexpiWarning(UserWarning):
    """Warns about information that could not be translated."""


class MappingError(ValueError):
    """Raised when a process step cannot be modeled by a unit class."""


def init_parameter_names(unit_class: type) -> tuple[str, ...]:
    """
    Return the names of the parameters that are specific to a unit class.

    Examples
    --------
    >>> import biosteam as bst
    >>> from biosteam.interop.dexpi.registry import init_parameter_names
    >>> init_parameter_names(bst.Mixer)
    ('rigorous', 'conserve_phases')

    """
    init = getattr(unit_class, '_init', None)
    if not inspect.isfunction(init):
        return ()
    kinds = (inspect.Parameter.POSITIONAL_OR_KEYWORD,
             inspect.Parameter.KEYWORD_ONLY)
    parameters = list(inspect.signature(init).parameters.values())[1:]
    return tuple(i.name for i in parameters if i.kind in kinds)


class UnitCodec:
    """
    Translate a family of unit operations to and from DEXPI process steps.

    The base class handles any unit operation generically: on export, the
    arguments of the unit's ``_init`` method are read back from attributes
    of the same name (BioSTEAM's convention) and stored in the BioSTEAM
    annotation of the process step. Subclasses refine this for parameters
    that are not plain data and add standard DEXPI attributes.

    """

    #: Names of ``_init`` parameters that are never stored.
    excluded_parameters: frozenset = frozenset()

    def dexpi_class(self, unit: Unit, default: str) -> str:
        """Return the DEXPI class representing `unit`."""
        return default

    def step_data(self, unit: Unit) -> dict[str, Any]:
        """Return the DEXPI data attributes describing `unit`."""
        return {}

    def step_quantities(self, unit: Unit) -> dict[str, tuple[float, str]]:
        """
        Return the DEXPI quantity attributes describing `unit` as pairs of
        a value in BioSTEAM base units and a quantity kind.
        """
        quantities = {}
        T = getattr(unit, 'T', None)
        if isinstance(T, Real) and not isinstance(T, bool):
            quantities['Temperature'] = (float(T), 'temperature')
        P = getattr(unit, 'P', None)
        if isinstance(P, Real) and not isinstance(P, bool):
            quantities['Pressure'] = (float(P), 'pressure')
        return quantities

    def annotate_step(self, unit: Unit, step: DexpiObject) -> None:
        """Add further attributes of `unit` to its process step."""

    def read_parameter(self, unit: Unit, name: str) -> Any:
        """
        Return the value of an ``_init`` parameter of `unit`.

        Raises
        ------
        AttributeError
            If the unit does not keep the parameter.

        """
        try:
            return getattr(unit, name)
        except AttributeError:
            return getattr(unit, '_' + name)

    def encode_parameters(self, unit: Unit) -> tuple[dict, list[str]]:
        """
        Return the parameters of `unit` as JSON-compatible values, and the
        names of the parameters that could not be stored.
        """
        parameters, omitted = {}, []
        for name in init_parameter_names(type(unit)):
            if name in self.excluded_parameters:
                continue
            try:
                parameters[name] = encode_value(
                    self.read_parameter(unit, name))
            except AttributeError:
                omitted.append(name)
            except ValueNotSerializable as error:
                omitted.append(name)
                warnings.warn(
                    f"parameter '{name}' of {unit!r} is not exported: {error}",
                    DexpiWarning, stacklevel=5,
                )
        return parameters, omitted

    def decode_parameters(self, parameters: Mapping[str, Any],
                          chemicals: tmo.Chemicals) -> dict:
        """Return keyword arguments for the unit from stored parameters."""
        return {i: decode_value(j, chemicals) for i, j in parameters.items()}

    def parameters_from_step(self, step: DexpiObject,
                             chemicals: tmo.Chemicals) -> dict:
        """
        Return keyword arguments for the unit from the standard DEXPI
        attributes of a process step without BioSTEAM annotation.

        Raises
        ------
        MappingError
            If the process step lacks information the unit class requires.

        """
        return {}

    def finalize(self, unit: Unit, step: Optional[DexpiObject] = None) -> None:
        """
        Complete a unit after it has been created by the importer from
        the process step `step`.
        """


@dataclass(frozen=True)
class Registration:
    """A unit class together with its DEXPI class and codec."""

    unit_class: type
    dexpi_class: str
    codec: UnitCodec


class UnitRegistry:
    """
    Registry of the mapping between BioSTEAM unit classes and DEXPI process
    step classes.

    Examples
    --------
    Lookups follow inheritance in both directions. A subclass of a
    registered unit class is exported like its parent:

    >>> import biosteam as bst
    >>> from biosteam.interop import dexpi
    >>> class MyPump(bst.Pump): pass
    >>> dexpi.default_registry.lookup_unit_class(MyPump).dexpi_class
    'Pumping'

    A DEXPI class without a registration of its own is imported like its
    nearest registered parent class:

    >>> dexpi.default_registry.lookup_dexpi_class(
    ...     'RotaryMixing').unit_class.__name__
    'Mixer'

    Register additional classes, optionally with a specialized codec:

    >>> registry = dexpi.UnitRegistry()
    >>> registry.register(MyPump, 'Pumping', import_default=True)
    >>> registry.lookup_dexpi_class('Pumping').unit_class.__name__
    'MyPump'

    """

    def __init__(self):
        self._by_unit_class: dict[type, Registration] = {}
        self._import_defaults: dict[str, Registration] = {}
        self._fallback = Registration(Unit, PROCESS_STEP, UnitCodec())

    def register(self, unit_class: type, dexpi_class: str,
                 codec: Optional[UnitCodec] = None, *,
                 import_default: bool = False,
                 import_default_for: tuple[str, ...] = ()) -> None:
        """
        Map a unit class to a DEXPI process step class.

        Parameters
        ----------
        unit_class :
            BioSTEAM unit class; subclasses are covered unless registered
            themselves.
        dexpi_class :
            Name of the DEXPI class written on export.
        codec :
            Translator of the unit's parameters. Defaults to the generic
            :class:`UnitCodec`.
        import_default :
            Whether steps of `dexpi_class` (and its subclasses) without
            BioSTEAM annotation are imported as `unit_class`.
        import_default_for :
            Further DEXPI classes to import as `unit_class` by default.

        """
        if not (isinstance(unit_class, type) and issubclass(unit_class, Unit)):
            raise TypeError("unit_class must be a subclass of Unit")
        registration = Registration(
            unit_class, dexpi_class, codec or UnitCodec())
        self._by_unit_class[unit_class] = registration
        if import_default:
            self._import_defaults[dexpi_class] = registration
        for name in import_default_for:
            self._import_defaults[name] = registration

    def lookup_unit_class(self, unit_class: type) -> Registration:
        """Return the registration that applies to a unit class."""
        for cls in unit_class.__mro__:
            registration = self._by_unit_class.get(cls)
            if registration is not None:
                return registration
        return self._fallback

    def lookup_dexpi_class(self, dexpi_class: str) -> Optional[Registration]:
        """Return the default registration for importing a DEXPI class."""
        for name in ancestors(dexpi_class):
            registration = self._import_defaults.get(name)
            if registration is not None:
                return registration
        return None

    def resolve_class_name(self, name: str) -> Optional[type]:
        """
        Return the unit class with a given name, or None.

        Only registered classes and the unit classes of the :mod:`biosteam`
        namespace are considered (BioSTEAM adds user-defined unit classes
        to its namespace when they are created). Names are never imported,
        so documents cannot trigger the execution of arbitrary modules.
        """
        for unit_class in reversed(self._by_unit_class):
            if unit_class.__name__ == name:
                return unit_class
        import biosteam as bst
        candidate = getattr(bst, name, None) if name.isidentifier() else None
        if isinstance(candidate, type) and issubclass(candidate, Unit):
            return candidate
        return None

    def copy(self) -> UnitRegistry:
        """Return an independent copy of the registry."""
        registry = UnitRegistry()
        registry._by_unit_class.update(self._by_unit_class)
        registry._import_defaults.update(self._import_defaults)
        return registry
