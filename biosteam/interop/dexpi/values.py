# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Encoding of unit operation parameters as JSON-compatible values.

Plain data (numbers, strings, lists and dictionaries) is stored as is.
Objects with richer meaning are wrapped in a dictionary with a ``__type__``
entry so that they can be restored: reactions, per-chemical data, and
non-finite numbers (which JSON cannot represent).

"""
from __future__ import annotations

import math
from numbers import Integral, Real
from typing import Any, Optional

import numpy as np
import thermosteam as tmo

from ..reactions import is_reaction, reaction_from_dict, reaction_to_dict

__all__ = ('ValueNotSerializable', 'decode_value', 'encode_value')

_TYPE = '__type__'


class ValueNotSerializable(TypeError):
    """Raised when a value cannot be stored in a document."""


def _is_chemical_data(value: Any) -> bool:
    return (
        hasattr(value, 'chemicals') and hasattr(value, 'data')
        and not hasattr(value, 'phases') and not is_reaction(value)
    )


def encode_value(value: Any) -> Any:
    """
    Return a JSON-compatible representation of a parameter value.

    Raises
    ------
    ValueNotSerializable
        If the value (or one of its items) has no representation.

    Examples
    --------
    >>> from biosteam.interop.dexpi.values import encode_value
    >>> encode_value((1, 2.5, 'a', None))
    [1, 2.5, 'a', None]
    >>> encode_value(float('inf'))
    {'__type__': 'float', 'value': 'inf'}

    """
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, Integral):
        return int(value)
    if isinstance(value, Real):
        value = float(value)
        if math.isfinite(value):
            return value
        return {_TYPE: 'float', 'value': repr(value)}
    if is_reaction(value):
        try:
            return {_TYPE: 'reaction', 'value': reaction_to_dict(value)}
        except NotImplementedError as error:
            raise ValueNotSerializable(str(error)) from None
    if _is_chemical_data(value):
        IDs = value.chemicals.IDs
        data = np.asarray(value.data, dtype=float).ravel()
        return {
            _TYPE: 'chemical data',
            'value': {IDs[i]: float(x) for i, x in enumerate(data) if x},
        }
    if isinstance(value, np.ndarray):
        return [encode_value(i) for i in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [encode_value(i) for i in value]
    if isinstance(value, dict):
        if not all(isinstance(i, str) for i in value) or _TYPE in value:
            raise ValueNotSerializable(
                "only dictionaries with string keys can be serialized"
            )
        return {i: encode_value(j) for i, j in value.items()}
    raise ValueNotSerializable(
        f"'{type(value).__name__}' objects cannot be serialized"
    )


def decode_value(value: Any, chemicals: Optional[tmo.Chemicals] = None) -> Any:
    """
    Return the parameter value of a representation created by
    :func:`encode_value`.

    Parameters
    ----------
    value :
        Encoded value.
    chemicals :
        Chemicals that reactions are defined for.

    Examples
    --------
    >>> from biosteam.interop.dexpi.values import decode_value
    >>> decode_value({'__type__': 'float', 'value': 'inf'})
    inf
    >>> decode_value({'__type__': 'chemical data', 'value': {'Water': 0.5}})
    {'Water': 0.5}

    """
    if isinstance(value, list):
        return [decode_value(i, chemicals) for i in value]
    if isinstance(value, dict):
        kind = value.get(_TYPE)
        if kind is None:
            return {i: decode_value(j, chemicals) for i, j in value.items()}
        content = value.get('value')
        if kind == 'float':
            return float(content)
        if kind == 'reaction':
            return reaction_from_dict(content, chemicals)
        if kind == 'chemical data':
            return dict(content)
        raise ValueError(f"unknown value type {kind!r}")
    return value
