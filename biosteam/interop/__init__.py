# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Interoperability with open data exchange standards.

* :mod:`biosteam.interop.dexpi` exchanges process flowsheets with DEXPI 2.0
  (Data Exchange in the Process Industry).
* :mod:`biosteam.interop.enzymeml` loads enzyme-catalyzed reactions and
  their kinetics from EnzymeML documents.
* :mod:`biosteam.interop.chemicals` resolves the chemical identifiers found
  in such documents to thermosteam chemicals.

The submodules are imported on first access, so importing BioSTEAM does not
depend on them.

"""
from __future__ import annotations

import importlib

__all__ = ('chemicals', 'dexpi', 'enzymeml', 'examples', 'reactions')


def __getattr__(name):
    if name in __all__:
        return importlib.import_module(f'{__name__}.{name}')
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def __dir__():
    return sorted(__all__)
