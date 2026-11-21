# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Exchange of process flowsheets with DEXPI 2.0.

DEXPI (Data Exchange in the Process Industry) is an open standard for
process engineering data. This module converts between BioSTEAM systems
and the DEXPI *Process* information model, which describes block flow and
process flow diagrams: process steps, their ports, the streams connecting
them, and the materials they carry.

.. autofunction:: biosteam.interop.dexpi.to_dexpi
.. autofunction:: biosteam.interop.dexpi.from_dexpi
.. autoclass:: biosteam.interop.dexpi.UnitRegistry
   :members:
.. autoclass:: biosteam.interop.dexpi.UnitCodec
   :members:
.. autoclass:: biosteam.interop.dexpi.GenericProcessStep
.. autoclass:: biosteam.interop.dexpi.DexpiExporter
   :members:
.. autoclass:: biosteam.interop.dexpi.DexpiImporter
   :members:

"""
from __future__ import annotations

import os
from datetime import datetime
from typing import Optional, Union

from . import (  # noqa: F401
    model, notes, quantities, validation, values, vocabulary, xmlio,
)
from .codecs import create_default_registry, default_registry
from .exporter import DexpiExporter
from .generic import GenericProcessStep
from .importer import DexpiImporter, DexpiImportError
from .model import DexpiDocument, DexpiObject
from .registry import DexpiWarning, MappingError, UnitCodec, UnitRegistry
from .validation import DexpiValidationError
from .xmlio import DexpiFormatError

__all__ = (
    'DexpiDocument',
    'DexpiExporter',
    'DexpiFormatError',
    'DexpiImportError',
    'DexpiImporter',
    'DexpiObject',
    'DexpiValidationError',
    'DexpiWarning',
    'GenericProcessStep',
    'MappingError',
    'UnitCodec',
    'UnitRegistry',
    'create_default_registry',
    'default_registry',
    'from_dexpi',
    'to_dexpi',
)


def to_dexpi(system, file: Union[str, os.PathLike, None] = None, *,
             registry: Optional[UnitRegistry] = None,
             timestamp: Optional[datetime] = None) -> Optional[str]:
    """
    Export a system as a DEXPI 2.0 process model.

    Parameters
    ----------
    system : System
        System to export.
    file :
        Path of the XML file to write. If not given, the XML text is
        returned instead.
    registry :
        Mapping between unit classes and DEXPI classes. Defaults to
        :data:`default_registry`.
    timestamp :
        Export time written to the document. Defaults to the current time.

    Returns
    -------
    str or None
        The document as XML text if no `file` is given.

    Notes
    -----
    Each unit operation becomes a process step of the DEXPI class it is
    registered for (a plain ``ProcessStep`` if it is not registered), each
    stream a DEXPI ``Stream`` with its temperature, pressure, flow and
    composition, and each chemical a material component. Feeds and
    products are terminated by ``Source`` and ``Sink`` steps. Parameters
    that DEXPI has no attributes for are stored in ``Note`` objects, from
    which :func:`from_dexpi` restores the system.

    Examples
    --------
    >>> import biosteam as bst
    >>> from biosteam.interop import dexpi
    >>> bst.main_flowsheet.set_flowsheet('to_dexpi_example')
    >>> bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> feed = bst.Stream('feed', Water=80, Ethanol=20)
    >>> P1 = bst.Pump('P1', ins=feed, P=2 * 101325)
    >>> S1 = bst.Splitter(
    ...     'S1', ins=P1-0, outs=('recycle', 'product'), split=0.2)
    >>> sys = bst.System.from_units('sys', [P1, S1])
    >>> sys.simulate()
    >>> xml = dexpi.to_dexpi(sys)
    >>> print(xml.splitlines()[1])
    <Model name="sys" uri="https://biosteam.readthedocs.io">

    The document can be loaded again as an equivalent system:

    >>> copy = dexpi.from_dexpi(xml, flowsheet='from_dexpi_example')
    >>> copy.simulate()
    >>> copy.flowsheet.stream.product.imol['Water']
    64.0

    """
    document = DexpiExporter(registry, timestamp).export(system)
    validation.validate(document)
    if file is None:
        return xmlio.dumps(document)
    xmlio.write(document, file)
    return None


def from_dexpi(source: Union[str, bytes, os.PathLike, DexpiDocument],
               ID: Optional[str] = None, *,
               chemicals=None,
               resolver=None,
               registry: Optional[UnitRegistry] = None,
               flowsheet=None,
               validate: bool = True):
    """
    Create a system from a DEXPI 2.0 process model.

    Parameters
    ----------
    source :
        Path of a DEXPI XML file, XML text, or a document object.
    ID :
        ID of the system. Defaults to the ID stored in the document.
    chemicals : Thermo or Iterable[Chemical], optional
        Chemicals (or a property package) to use for the material
        components of the document, matched by ID, CAS number or InChIKey.
    resolver : ResolverChain, optional
        Resolver for material components that `chemicals` does not cover.
        Defaults to :func:`biosteam.interop.chemicals.default_resolver`,
        which searches thermosteam's databases and creates placeholders
        with default properties for unknown components.
    registry :
        Mapping between unit classes and DEXPI classes. Defaults to
        :data:`default_registry`.
    flowsheet : Flowsheet or str, optional
        Flowsheet (or the ID of a new one) in which the units and streams
        are registered. Defaults to the main flowsheet.
    validate :
        Whether to check the structure of the document first.

    Returns
    -------
    System
        The system, ready to be simulated.

    Notes
    -----
    If neither `chemicals` nor `resolver` is given and the active property
    package (``settings.get_thermo()``) defines exactly the components of
    the document, it is used as is. If no property package is active, the
    one created for the document becomes the active one.

    Process steps written by BioSTEAM are restored with their original
    class and parameters. Other process steps are created from the unit
    class registered for their DEXPI class, or as
    :class:`GenericProcessStep` if there is none, so that the topology is
    always complete.

    Raises
    ------
    DexpiFormatError
        If `source` is not well-formed DEXPI XML.
    DexpiValidationError
        If the process model is structurally invalid.
    DexpiImportError
        If the document cannot be interpreted as a system.

    """
    document = (
        source if isinstance(source, DexpiDocument) else xmlio.read(source)
    )
    if validate:
        validation.validate(document)
    importer = DexpiImporter(registry, chemicals, resolver, flowsheet)
    return importer.load(document, ID)
