# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Conversion of thermosteam reactions to and from plain dictionaries.

Reactions are part of a unit operation's specification but have no
counterpart in flowsheet exchange formats. Turning them into dictionaries
of JSON-compatible values lets them travel as annotations.

Kinetic reactions carry a `rate` callable. A rate can be serialized if it
implements the :class:`SerializableRate` protocol; its type name must be
registered with :func:`register_rate_type` to be restored.

"""
from __future__ import annotations

import importlib
import math
import warnings
from typing import (
    Any, Callable, Iterator, Mapping, Optional, Protocol, runtime_checkable,
)

import thermosteam as tmo

__all__ = (
    'ReactionSerializationWarning',
    'SerializableRate',
    'is_reaction',
    'iter_rates',
    'reaction_from_dict',
    'reaction_to_dict',
    'register_rate_type',
)

_REACTION_SETS = {
    'ParallelReaction': tmo.ParallelReaction,
    'SeriesReaction': tmo.SeriesReaction,
}

#: Loaders of rate objects by type name. Values are either callables or
#: 'module:attribute' paths of loaders within BioSTEAM, imported on demand.
_rate_loaders: dict[str, Any] = {
    'SteadyStateCSTRRate':
        'biosteam.interop.enzymeml.kinetics:SteadyStateCSTRRate.from_dict',
}


class ReactionSerializationWarning(UserWarning):
    """Warns that part of a reaction could not be serialized or restored."""


@runtime_checkable
class SerializableRate(Protocol):
    """Protocol of reaction rates that can be stored in a document."""

    def to_dict(self) -> dict:
        """Return the rate as a dictionary with a 'type' entry."""


def register_rate_type(
        name: str, loader: Callable[[Mapping, tmo.Chemicals], Callable]
) -> None:
    """
    Register a loader that restores rates of a given type.

    Parameters
    ----------
    name :
        Value of the 'type' entry written by the rate's `to_dict` method.
    loader :
        Function that accepts the dictionary and the chemicals and returns
        the rate.

    """
    _rate_loaders[name] = loader


def _get_rate_loader(name: str) -> Optional[Callable]:
    loader = _rate_loaders.get(name)
    if isinstance(loader, str):
        module_name, _, attribute = loader.partition(':')
        loader = importlib.import_module(module_name)
        for part in attribute.split('.'):
            loader = getattr(loader, part)
        _rate_loaders[name] = loader
    return loader


def is_reaction(obj: Any) -> bool:
    """Return whether `obj` is a thermosteam reaction of any kind."""
    return isinstance(
        obj, (tmo.Reaction, tmo.reaction.ReactionSet, tmo.ReactionSystem))


def iter_rates(reaction: Any) -> Iterator[Callable]:
    """Iterate over the kinetic rates of all reactions within `reaction`."""
    if reaction is None:
        return
    if isinstance(reaction, tmo.ReactionSystem):
        for i in reaction._reactions:
            yield from iter_rates(i)
    elif isinstance(reaction, tmo.reaction.ReactionSet):
        for i in reaction:
            yield from iter_rates(i)
    elif isinstance(reaction, tmo.Reaction):
        if reaction.rate is not None:
            yield reaction.rate


def _stoichiometry(reaction: tmo.Reaction) -> dict[str, float]:
    IDs = reaction.chemicals.IDs
    stoichiometry = reaction.stoichiometry
    return {
        IDs[i]: float(stoichiometry[i])
        for i in range(len(IDs)) if stoichiometry[i]
    }


def _single_to_dict(reaction: tmo.Reaction) -> dict:
    if reaction._phases:
        raise NotImplementedError(
            "reactions with phase-specific stoichiometry cannot be serialized"
        )
    data = {
        'type': 'Reaction',
        'stoichiometry': _stoichiometry(reaction),
        'reactant': reaction.reactant,
        'basis': reaction.basis,
    }
    rate = reaction.rate
    if rate is None:
        data['X'] = float(reaction.X)
    elif isinstance(rate, SerializableRate):
        data['rate'] = rate.to_dict()
    else:
        warnings.warn(
            f"the rate of reaction '{reaction.reactant}' is a "
            f"'{type(rate).__name__}' object, which cannot be serialized; "
            "the reaction is stored without kinetics",
            ReactionSerializationWarning, stacklevel=4,
        )
        data['rate'] = None
    return data


def reaction_to_dict(reaction: Any) -> dict:
    """
    Return a reaction, reaction set or reaction system as a dictionary.

    Examples
    --------
    >>> import thermosteam as tmo
    >>> from biosteam.interop.reactions import reaction_to_dict
    >>> chemicals = tmo.Chemicals(['H2', 'O2', 'Water'], cache=True)
    >>> chemicals.compile()
    >>> reaction = tmo.Reaction('2H2 + O2 -> 2Water', 'H2', 0.9, chemicals)
    >>> reaction_to_dict(reaction)
    {'type': 'Reaction',
     'stoichiometry': {'H2': -1.0, 'O2': -0.5, 'Water': 1.0},
     'reactant': 'H2', 'basis': 'mol', 'X': 0.9}

    """
    if isinstance(reaction, tmo.ReactionSystem):
        return {
            'type': 'ReactionSystem',
            'reactions': [reaction_to_dict(i) for i in reaction._reactions],
        }
    if isinstance(reaction, tmo.reaction.ReactionSet):
        for name, cls in _REACTION_SETS.items():
            if isinstance(reaction, cls):
                break
        else:
            raise NotImplementedError(
                f"'{type(reaction).__name__}' objects cannot be serialized"
            )
        return {
            'type': name,
            'reactions': [_single_to_dict(i) for i in reaction],
        }
    if isinstance(reaction, tmo.Reaction):
        return _single_to_dict(reaction)
    raise TypeError(f"'{type(reaction).__name__}' object is not a reaction")


def _single_from_dict(data: Mapping, chemicals: tmo.Chemicals) -> tmo.Reaction:
    rate = None
    X = data.get('X')
    if 'rate' in data:
        rate_data = data['rate']
        loader = _get_rate_loader(rate_data.get('type')) if rate_data else None
        if loader is None:
            warnings.warn(
                f"the kinetics of reaction '{data.get('reactant')}' could "
                "not be restored; assuming no conversion until a rate is set",
                ReactionSerializationWarning, stacklevel=4,
            )
            X = 0.
        else:
            rate = loader(rate_data, chemicals)
    if X is not None and isinstance(X, float) and math.isnan(X):
        X = 0.
    return tmo.Reaction(
        dict(data['stoichiometry']), data['reactant'], X,
        chemicals=chemicals, basis=data.get('basis', 'mol'), rate=rate,
    )


def reaction_from_dict(data: Mapping, chemicals: tmo.Chemicals) -> Any:
    """
    Return the reaction described by a dictionary.

    Parameters
    ----------
    data :
        Dictionary created by :func:`reaction_to_dict`.
    chemicals :
        Chemicals the reaction is defined for.

    Examples
    --------
    >>> import thermosteam as tmo
    >>> from biosteam.interop.reactions import reaction_from_dict
    >>> chemicals = tmo.Chemicals(['H2', 'O2', 'Water'], cache=True)
    >>> chemicals.compile()
    >>> reaction = reaction_from_dict(
    ...     {'type': 'Reaction', 'reactant': 'H2', 'X': 0.9,
    ...      'stoichiometry': {'H2': -1.0, 'O2': -0.5, 'Water': 1.0}},
    ...     chemicals)
    >>> reaction.show()
    Reaction (by mol):
    stoichiometry             reactant    X[%]
    H2 + 0.5 O2 -> Water      H2         90.00

    """
    kind = data.get('type', 'Reaction')
    if kind == 'ReactionSystem':
        return tmo.ReactionSystem(
            *[reaction_from_dict(i, chemicals) for i in data['reactions']]
        )
    if kind in _REACTION_SETS:
        return _REACTION_SETS[kind](
            [_single_from_dict(i, chemicals) for i in data['reactions']]
        )
    if kind == 'Reaction':
        return _single_from_dict(data, chemicals)
    raise ValueError(f"unknown reaction type {kind!r}")
