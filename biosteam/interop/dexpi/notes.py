# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
BioSTEAM annotations in DEXPI documents.

DEXPI has no attributes for simulator-specific data such as split
fractions or reactions. BioSTEAM stores them as JSON in standard DEXPI
``Note`` objects classified as ``'BioSTEAM'``, which any DEXPI tool can
read, display and write back unchanged.

"""
from __future__ import annotations

import json
from typing import Optional

from .model import NOTE, DexpiObject, multi_language_string, read_text

__all__ = ('NOTE_CLASSIFICATION', 'SCHEMA_VERSION', 'create_note', 'read_note')

#: Classification of notes written by BioSTEAM.
NOTE_CLASSIFICATION = 'BioSTEAM'

#: Version of the annotation layout.
SCHEMA_VERSION = 1


def create_note(kind: str, payload: dict) -> DexpiObject:
    """
    Return a note holding a BioSTEAM annotation.

    Parameters
    ----------
    kind :
        What the annotated object is: 'system', 'unit', 'stream' or
        'chemical'.
    payload :
        JSON-compatible content of the annotation.

    """
    content = {'schema': SCHEMA_VERSION, 'kind': kind, **payload}
    text = json.dumps(content, sort_keys=True, allow_nan=False)
    return DexpiObject(NOTE, data={
        'NoteClassification': NOTE_CLASSIFICATION,
        'NoteText': multi_language_string(text),
    })


def read_note(obj: DexpiObject) -> Optional[dict]:
    """
    Return the BioSTEAM annotation of an object, or None if it has none.

    Notes that are not classified as BioSTEAM annotations, or that do not
    hold a JSON object, are ignored.

    Examples
    --------
    >>> from biosteam.interop.dexpi import model, notes
    >>> step = model.DexpiObject(model.process_type('Pumping'))
    >>> step.refer('ReferencedNotes',
    ...            notes.create_note('unit', {'class': 'Pump'}))
    >>> notes.read_note(step)
    {'class': 'Pump', 'kind': 'unit', 'schema': 1}

    """
    for note in obj.targets('ReferencedNotes'):
        if note.get('NoteClassification') != NOTE_CLASSIFICATION:
            continue
        text = read_text(note.get('NoteText'))
        if not text:
            continue
        try:
            content = json.loads(text)
        except ValueError:
            continue
        if isinstance(content, dict):
            return content
    return None
