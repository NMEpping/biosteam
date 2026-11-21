# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Reader and writer for the DEXPI 2.0 XML exchange format.

The format has a small, regular grammar::

    <Model name=".." uri="..">
      <Import prefix="Core" source=".../Core.xml"/>
      <Object id=".." type="Core/EngineeringModel">
        <Components property="..."> <Object .../> ... </Components>
        <References property="..." objects="#id1 #id2"/>
        <Data property="..."> <String>..</String> </Data>
      </Object>
    </Model>

Only the standard library is used, so reading and writing DEXPI documents
adds no dependency to BioSTEAM.

"""
from __future__ import annotations

import math
import os
import re
import xml.etree.ElementTree as ET
from datetime import datetime
from typing import Any, Union

from .model import (
    DEXPI_VERSION,
    MODEL_URI,
    PACKAGES,
    UNDEFINED,
    Aggregate,
    DexpiDocument,
    DexpiObject,
    Enumeration,
)

__all__ = ('DexpiFormatError', 'dumps', 'loads', 'read', 'write')

_INVALID_XML_CHARACTERS = re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f]')
_FORBIDDEN_DECLARATIONS = re.compile(rb'<!\s*(DOCTYPE|ENTITY)', re.IGNORECASE)


class DexpiFormatError(ValueError):
    """Raised when a document is not well-formed DEXPI XML."""


# %% Writing

def _double_text(value: float) -> str:
    """Return the ``xsd:double`` literal of a float."""
    if math.isfinite(value):
        return repr(float(value))
    if math.isnan(value):
        return 'NaN'
    return 'INF' if value > 0 else '-INF'


def _write_value(parent: ET.Element, value: Any, location: str) -> None:
    """Append the XML element of a data value to `parent`."""
    if value is UNDEFINED:
        ET.SubElement(parent, 'Undefined')
    elif isinstance(value, bool):
        ET.SubElement(parent, 'Boolean').text = 'true' if value else 'false'
    elif isinstance(value, int):
        ET.SubElement(parent, 'Integer').text = str(value)
    elif isinstance(value, float):
        ET.SubElement(parent, 'Double').text = _double_text(value)
    elif isinstance(value, str):
        if _INVALID_XML_CHARACTERS.search(value):
            raise DexpiFormatError(
                f"{location} contains a control character that XML cannot "
                "represent"
            )
        ET.SubElement(parent, 'String').text = value
    elif isinstance(value, datetime):
        ET.SubElement(parent, 'DateTime').text = value.isoformat()
    elif isinstance(value, Enumeration):
        ET.SubElement(parent, 'DataReference', data=value.path)
    elif isinstance(value, Aggregate):
        element = ET.SubElement(parent, 'AggregatedDataValue', type=value.type)
        _write_data(element, value.data, location)
    else:
        raise DexpiFormatError(
            f"{location} has a value of type '{type(value).__name__}', "
            "which is not a DEXPI data value"
        )


def _write_data(parent: ET.Element, data: dict[str, Any],
                location: str) -> None:
    """Append the ``Data`` elements of all data attributes to `parent`."""
    for name, value in data.items():
        if value is None:
            continue  # An absent attribute and None are equivalent.
        element = ET.SubElement(parent, 'Data', property=name)
        values = value if isinstance(value, list) else [value]
        for item in values:
            _write_value(element, item, f"attribute '{name}' of {location}")


def _assign_ids(root: DexpiObject) -> None:
    """Give every object of the tree a unique identifier."""
    taken = {i.id for i in root.walk() if i.id}
    counter = 0
    for obj in root.walk():
        if obj.id:
            continue
        while True:
            counter += 1
            candidate = f'id_{counter}'
            if candidate not in taken:
                break
        obj.id = candidate
        taken.add(candidate)


def _write_object(parent: ET.Element, obj: DexpiObject) -> None:
    """Append the XML element of an object to `parent`."""
    element = ET.SubElement(parent, 'Object', id=obj.id, type=obj.type)
    for name, children in obj.components.items():
        if not children:
            continue
        components = ET.SubElement(element, 'Components', property=name)
        for child in children:
            _write_object(components, child)
    for name, targets in obj.references.items():
        if not targets:
            continue
        missing = [i for i in targets if not i.id]
        if missing:
            raise DexpiFormatError(
                f"reference '{name}' of {obj!r} points to {missing[0]!r}, "
                "which is not part of the document"
            )
        ET.SubElement(
            element, 'References',
            objects=' '.join(f'#{i.id}' for i in targets),
            property=name,
        )
    _write_data(element, obj.data, repr(obj))


def to_element(document: DexpiDocument) -> ET.Element:
    """Return the ``Model`` XML element of a document."""
    _assign_ids(document.root)
    model = ET.Element('Model', name=document.name, uri=document.uri)
    for package in PACKAGES:
        ET.SubElement(
            model, 'Import', prefix=package,
            source=MODEL_URI.format(version=DEXPI_VERSION, package=package),
        )
    _write_object(model, document.root)
    return model


def dumps(document: DexpiDocument) -> str:
    """
    Return a document as DEXPI XML text.

    Examples
    --------
    >>> from biosteam.interop.dexpi import model, xmlio
    >>> root = model.DexpiObject(model.ENGINEERING_MODEL)
    >>> root.data['OriginatingSystemName'] = 'BioSTEAM'
    >>> print(xmlio.dumps(model.DexpiDocument(root, name='demo')))
    <?xml version='1.0' encoding='utf-8'?>
    <Model name="demo" uri="https://biosteam.readthedocs.io">
      <Import prefix="Core"
              source="https://data.dexpi.org/models/2.0.1/Core.xml" />
      <Import prefix="Process"
              source="https://data.dexpi.org/models/2.0.1/Process.xml" />
      <Object id="id_1" type="Core/EngineeringModel">
        <Data property="OriginatingSystemName">
          <String>BioSTEAM</String>
        </Data>
      </Object>
    </Model>

    """
    element = to_element(document)
    ET.indent(element, space='  ')
    text = ET.tostring(element, encoding='unicode')
    return f"<?xml version='1.0' encoding='utf-8'?>\n{text}\n"


def write(document: DexpiDocument, file: Union[str, os.PathLike]) -> None:
    """Write a document to a DEXPI XML file."""
    with open(file, 'w', encoding='utf-8', newline='\n') as f:
        f.write(dumps(document))


# %% Reading

def _read_value(element: ET.Element) -> Any:
    """Return the data value held by a value element."""
    tag = element.tag
    text = element.text or ''
    if tag == 'String':
        return text
    if tag == 'Double':
        try:
            return float(text.strip().replace('INF', 'inf'))
        except ValueError:
            raise DexpiFormatError(f"invalid Double value {text!r}") from None
    if tag == 'Integer':
        try:
            return int(text.strip())
        except ValueError:
            raise DexpiFormatError(f"invalid Integer value {text!r}") from None
    if tag == 'Boolean':
        return text.strip().lower() in ('true', '1')
    if tag == 'DateTime':
        try:
            return datetime.fromisoformat(text.strip())
        except ValueError:
            return text.strip()
    if tag == 'DataReference':
        return Enumeration(element.get('data', ''))
    if tag == 'AggregatedDataValue':
        return Aggregate(element.get('type', ''), _read_data(element))
    if tag == 'Undefined':
        return UNDEFINED
    if tag == 'Null':
        return None
    # Unknown value kinds are kept as text so that no information is lost.
    return text


def _read_data(parent: ET.Element) -> dict[str, Any]:
    """Return the data attributes declared directly under `parent`."""
    data = {}
    for element in parent.findall('Data'):
        name = element.get('property')
        if not name:
            raise DexpiFormatError("Data element without a 'property'")
        values = [_read_value(i) for i in element]
        if name in data:
            previous = data[name]
            if not isinstance(previous, list):
                previous = [previous]
            data[name] = previous + values
        elif len(values) == 1:
            data[name] = values[0]
        else:
            data[name] = values
    return data


def _read_object(element: ET.Element, index: dict, pending: list,
                 prefixes: dict[str, str]) -> DexpiObject:
    """Return the object of an ``Object`` element, recursively."""
    type_ = element.get('type')
    if not type_:
        raise DexpiFormatError("Object element without a 'type'")
    prefix, separator, name = type_.partition('/')
    if separator and prefix in prefixes:
        type_ = f'{prefixes[prefix]}/{name}'
    obj = DexpiObject(type_, element.get('id'), _read_data(element))
    if obj.id:
        if obj.id in index:
            raise DexpiFormatError(f"duplicate object id '{obj.id}'")
        index[obj.id] = obj
    for components in element.findall('Components'):
        name = components.get('property')
        if not name:
            raise DexpiFormatError("Components element without a 'property'")
        children = [
            _read_object(i, index, pending, prefixes)
            for i in components.findall('Object')
        ]
        obj.components.setdefault(name, []).extend(children)
    for references in element.findall('References'):
        name = references.get('property')
        if not name:
            raise DexpiFormatError("References element without a 'property'")
        ids = [i.lstrip('#') for i in references.get('objects', '').split()]
        pending.append((obj, name, ids))
    return obj


def _canonical_prefixes(model: ET.Element) -> dict[str, str]:
    """Map the prefixes declared by a document to DEXPI package names."""
    prefixes = {}
    for element in model.findall('Import'):
        prefix = element.get('prefix')
        source = element.get('source', '')
        package = source.rstrip('/').rsplit('/', 1)[-1]
        if package.endswith('.xml'):
            package = package[:-4]
        if prefix and package:
            prefixes[prefix] = package
    return prefixes


def loads(text: Union[str, bytes]) -> DexpiDocument:
    """
    Return the document described by DEXPI XML text.

    Raises
    ------
    DexpiFormatError
        If the text is not well-formed DEXPI XML, declares a document type
        or entities, or holds a reference to an object it does not define.

    Examples
    --------
    >>> from biosteam.interop.dexpi import xmlio
    >>> document = xmlio.loads('''
    ... <Model name="demo" uri="urn:demo">
    ...   <Object id="a" type="Core/EngineeringModel">
    ...     <Data property="OriginatingSystemName"><String>X</String></Data>
    ...   </Object>
    ... </Model>''')
    >>> document.root.get('OriginatingSystemName')
    'X'

    """
    raw = text.encode('utf-8') if isinstance(text, str) else bytes(text)
    if _FORBIDDEN_DECLARATIONS.search(raw):
        # DEXPI documents never declare a DTD; refusing them rules out
        # entity-expansion and external-entity attacks.
        raise DexpiFormatError(
            "document type and entity declarations are not allowed in "
            "DEXPI documents"
        )
    try:
        model = ET.fromstring(raw)
    except ET.ParseError as error:
        raise DexpiFormatError(
            f"document is not well-formed XML: {error}") from None
    if model.tag != 'Model':
        raise DexpiFormatError(
            f"expected a 'Model' root element, found '{model.tag}'; "
            "only the DEXPI 2.0 XML format is supported"
        )
    roots = model.findall('Object')
    if len(roots) != 1:
        raise DexpiFormatError(
            f"expected exactly one root object, found {len(roots)}"
        )
    index: dict[str, DexpiObject] = {}
    pending: list = []
    root = _read_object(roots[0], index, pending, _canonical_prefixes(model))
    for obj, name, ids in pending:
        try:
            targets = [index[i] for i in ids]
        except KeyError as error:
            raise DexpiFormatError(
                f"reference '{name}' of {obj!r} points to undefined object "
                f"{error.args[0]!r}"
            ) from None
        obj.references.setdefault(name, []).extend(targets)
    return DexpiDocument(
        root,
        name=model.get('name', ''),
        uri=model.get('uri', ''),
    )


def read(source: Union[str, bytes, os.PathLike]) -> DexpiDocument:
    """
    Return the document of a DEXPI XML file or XML text.

    Parameters
    ----------
    source :
        Path of a file, or the XML text itself (as bytes, or as a string
        starting with ``<``).

    """
    if isinstance(source, bytes):
        return loads(source)
    if isinstance(source, str) and source.lstrip().startswith('<'):
        return loads(source)
    with open(source, 'rb') as f:
        return loads(f.read())
