# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
In-memory representation of a DEXPI 2.0 document.

A DEXPI document is a tree of objects. Each object has a type, an optional
identifier, and three kinds of attributes:

* *data* attributes hold values (strings, numbers, enumeration literals or
  aggregated values such as a physical quantity);
* *composition* attributes hold child objects owned by the object;
* *reference* attributes point to objects owned elsewhere in the tree.

:class:`DexpiObject` mirrors that structure one to one, so that reading and
writing a document is lossless regardless of how much of the DEXPI
vocabulary BioSTEAM understands. The names used in this module (types,
attributes and enumeration literals) are those of the DEXPI Specification
2.0, (c) DEXPI e.V., licensed under CC BY 4.0.

"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterator, Optional

__all__ = (
    'Aggregate',
    'DexpiDocument',
    'DexpiObject',
    'Enumeration',
    'UNDEFINED',
    'core_type',
    'multi_language_string',
    'physical_quantity',
    'process_type',
    'qualified_value',
    'read_quantity',
    'read_text',
)

#: Version of the DEXPI specification written by the exporter.
DEXPI_VERSION = '2.0.1'

#: Location of the DEXPI model definitions imported by a document.
MODEL_URI = 'https://data.dexpi.org/models/{version}/{package}.xml'

#: Packages declared by documents written by BioSTEAM.
PACKAGES = ('Core', 'Process')

ENGINEERING_MODEL = 'Core/EngineeringModel'
PROCESS_MODEL = 'Process/ProcessModel'
QUALIFIED_VALUE = 'Core/QualifiedValue'
NOTE = 'Core/Note'
PERSISTENT_IDENTIFIER = 'Core/PersistentIdentifier'
PHYSICAL_QUANTITY = 'Core/PhysicalQuantities.PhysicalQuantity'
MULTI_LANGUAGE_STRING = 'Core/DataTypes.MultiLanguageString'
SINGLE_LANGUAGE_STRING = 'Core/DataTypes.SingleLanguageString'


class _Undefined:
    """Type of the :data:`UNDEFINED` singleton."""

    __slots__ = ()
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __repr__(self):
        return 'UNDEFINED'

    def __bool__(self):
        return False


#: Value of a data attribute that is explicitly declared as undefined.
UNDEFINED = _Undefined()


@dataclass(frozen=True)
class Enumeration:
    """
    A literal of a DEXPI enumeration.

    Parameters
    ----------
    path :
        Qualified name of the literal, e.g.
        ``'Process/Enumerations.PortDirection.Inlet'``.

    Examples
    --------
    >>> from biosteam.interop.dexpi.model import Enumeration
    >>> direction = Enumeration('Process/Enumerations.PortDirection.Inlet')
    >>> direction.literal
    'Inlet'
    >>> direction.enumeration
    'PortDirection'

    """

    path: str

    @property
    def literal(self) -> str:
        """Name of the literal without its enumeration."""
        return self.path.rsplit('.', 1)[-1]

    @property
    def enumeration(self) -> str:
        """Name of the enumeration the literal belongs to."""
        parts = self.path.replace('/', '.').split('.')
        return parts[-2] if len(parts) > 1 else ''


@dataclass
class Aggregate:
    """
    An aggregated data value: a typed bundle of data attributes that has
    no identity of its own (e.g. a physical quantity).

    Parameters
    ----------
    type :
        Qualified name of the data type.
    data :
        Data attributes of the value by name.

    """

    type: str
    data: dict[str, Any] = field(default_factory=dict)


@dataclass(eq=False)
class DexpiObject:
    """
    A DEXPI object with data, composition and reference attributes.

    Parameters
    ----------
    type :
        Qualified type name, e.g. ``'Process/Process.Stream'``.
    id :
        Identifier used for references within the document. Assigned by
        the writer when missing.
    data :
        Data attributes by name. A list holds a multi-valued attribute.
    components :
        Owned child objects by composition attribute name.
    references :
        Referenced objects by reference attribute name.

    Examples
    --------
    >>> from biosteam.interop.dexpi.model import DexpiObject, process_type
    >>> stream = DexpiObject(process_type('Stream'), data={'Identifier': 's1'})
    >>> stream.class_name
    'Stream'
    >>> stream.get('Identifier')
    's1'
    >>> stream.get('Label', 'unnamed')
    'unnamed'

    """

    type: str
    id: Optional[str] = None
    data: dict[str, Any] = field(default_factory=dict)
    components: dict[str, list[DexpiObject]] = field(default_factory=dict)
    references: dict[str, list[DexpiObject]] = field(default_factory=dict)

    @property
    def class_name(self) -> str:
        """Type name without package qualification."""
        return self.type.rsplit('/', 1)[-1].rsplit('.', 1)[-1]

    def get(self, name: str, default: Any = None) -> Any:
        """Return a data attribute, or `default` if missing or undefined."""
        value = self.data.get(name, default)
        return default if value is UNDEFINED else value

    def children(self, name: str) -> list[DexpiObject]:
        """Return all objects of a composition attribute."""
        return self.components.get(name, [])

    def child(self, name: str) -> Optional[DexpiObject]:
        """Return the first object of a composition attribute, if any."""
        children = self.components.get(name)
        return children[0] if children else None

    def targets(self, name: str) -> list[DexpiObject]:
        """Return all objects of a reference attribute."""
        return self.references.get(name, [])

    def target(self, name: str) -> Optional[DexpiObject]:
        """Return the first object of a reference attribute, if any."""
        targets = self.references.get(name)
        return targets[0] if targets else None

    def add(self, name: str, *objects: DexpiObject) -> DexpiObject:
        """Append objects to a composition attribute and return the last."""
        self.components.setdefault(name, []).extend(objects)
        return objects[-1]

    def refer(self, name: str, *objects: DexpiObject) -> None:
        """Append objects to a reference attribute."""
        self.references.setdefault(name, []).extend(objects)

    def walk(self) -> Iterator[DexpiObject]:
        """Iterate over this object and all objects it owns, depth first."""
        yield self
        for children in self.components.values():
            for child in children:
                yield from child.walk()

    def __repr__(self):
        identifier = self.data.get('Identifier', self.id)
        return f"<{self.class_name}: {identifier}>"


@dataclass
class DexpiDocument:
    """
    A DEXPI document: a root object plus the document metadata.

    Parameters
    ----------
    root :
        Root object, normally a ``Core/EngineeringModel``.
    name :
        Name of the model.
    uri :
        URI identifying the model.

    """

    root: DexpiObject
    name: str = 'BioSTEAM'
    uri: str = 'https://biosteam.readthedocs.io'

    @property
    def conceptual_model(self) -> Optional[DexpiObject]:
        """The conceptual model of the document, if any."""
        if self.root.type == ENGINEERING_MODEL:
            return self.root.child('ConceptualModel')
        return self.root

    @property
    def process_model(self) -> Optional[DexpiObject]:
        """The process model of the document, or None if it has none."""
        model = self.conceptual_model
        if model is not None and model.class_name == 'ProcessModel':
            return model
        return None

    def objects(self) -> Iterator[DexpiObject]:
        """Iterate over all objects of the document."""
        return self.root.walk()


def process_type(name: str) -> str:
    """Return the qualified type of a class of the DEXPI process package."""
    return f'Process/Process.{name}'


def core_type(name: str) -> str:
    """Return the qualified type of a class of the DEXPI core package."""
    return f'Core/{name}'


def process_enumeration(enumeration: str, literal: str) -> Enumeration:
    """Return a literal of an enumeration of the DEXPI process package."""
    return Enumeration(f'Process/Enumerations.{enumeration}.{literal}')


def multi_language_string(text: str,
                          language: Optional[str] = None) -> Aggregate:
    """Return a multi-language string holding a single value."""
    data = {'Value': text}
    if language:
        data = {'Language': language, 'Value': text}
    return Aggregate(
        MULTI_LANGUAGE_STRING,
        {'SingleLanguageStrings': [Aggregate(SINGLE_LANGUAGE_STRING, data)]},
    )


def read_text(value: Any, default: Optional[str] = None) -> Optional[str]:
    """
    Return the text of a string or multi-language string data value.

    Examples
    --------
    >>> from biosteam.interop.dexpi.model import (
    ...     multi_language_string, read_text)
    >>> read_text(multi_language_string('reactor'))
    'reactor'
    >>> read_text('reactor')
    'reactor'

    """
    if isinstance(value, str):
        return value
    if isinstance(value, Aggregate):
        strings = value.data.get('SingleLanguageStrings', [])
        if not isinstance(strings, list):
            strings = [strings]
        for string in strings:
            if isinstance(string, Aggregate):
                text = string.data.get('Value')
                if isinstance(text, str):
                    return text
    return default


def physical_quantity(value: float, unit: str, literal: str) -> Aggregate:
    """
    Return a physical quantity data value.

    Parameters
    ----------
    value :
        Numerical value.
    unit :
        Name of the unit enumeration, e.g. ``'TemperatureUnit'``.
    literal :
        Name of the unit literal, e.g. ``'Kelvin'``.

    """
    return Aggregate(
        PHYSICAL_QUANTITY,
        {'Unit': Enumeration(f'Core/PhysicalQuantities.{unit}.{literal}'),
         'Value': float(value)},
    )


def qualified_value(value: Any) -> DexpiObject:
    """Return a qualified value object wrapping a data value."""
    return DexpiObject(QUALIFIED_VALUE, data={'Value': value})


def read_quantity(value: Any) -> Optional[tuple[float, Optional[Enumeration]]]:
    """
    Return the number and unit of a qualified value, physical quantity or
    plain number, or None if `value` holds no number.

    """
    if isinstance(value, DexpiObject):
        value = value.get('Value')
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value), None
    if isinstance(value, Aggregate):
        number = value.data.get('Value')
        unit = value.data.get('Unit')
        if isinstance(number, (int, float)) and not isinstance(number, bool):
            return (
                float(number), unit if isinstance(unit, Enumeration) else None
            )
    return None
