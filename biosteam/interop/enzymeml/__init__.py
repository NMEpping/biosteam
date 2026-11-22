# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Enzyme kinetics from EnzymeML documents.

EnzymeML is an open exchange format for enzyme-catalyzed reactions. This
module reads the reaction model of an EnzymeML document (species,
stoichiometry, rate law and kinetic parameters) and turns it into
thermosteam chemicals and kinetic reactions for stirred tank reactors.

EnzymeML JSON files are read without additional dependencies. OMEX
archives are read with the optional `pyenzyme` package.

.. autofunction:: biosteam.interop.enzymeml.read
.. autoclass:: biosteam.interop.enzymeml.EnzymeMLDocument
   :members:
.. autoclass:: biosteam.interop.enzymeml.SteadyStateCSTRRate
   :members:
.. autoclass:: biosteam.interop.enzymeml.PlugFlowRate
   :members:
.. autoclass:: biosteam.interop.enzymeml.RateLaw
   :members:

"""
from . import document, kinetics, ratelaw, units  # noqa: F401
from .document import (
    EnzymeMLDocument,
    EnzymeMLError,
    Parameter,
    Reaction,
    ReactionElement,
    Species,
    read,
)
from .kinetics import (
    PlugFlowRate,
    SteadyStateCSTRRate,
    create_chemicals,
    create_reaction,
)
from .ratelaw import RateLaw, RateLawError
from .units import UnitError, conversion_factor

__all__ = (
    'EnzymeMLDocument',
    'EnzymeMLError',
    'Parameter',
    'PlugFlowRate',
    'RateLaw',
    'RateLawError',
    'Reaction',
    'ReactionElement',
    'Species',
    'SteadyStateCSTRRate',
    'UnitError',
    'conversion_factor',
    'create_chemicals',
    'create_reaction',
    'read',
)
