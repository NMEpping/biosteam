# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Reference case for the interoperability modules: the regeneration of NADH
by formate dehydrogenase (FDH) in a continuous stirred tank reactor.

Formate serves as a sacrificial cosubstrate that reduces NAD+ to NADH:

.. math::

   \\text{Formate} + \\text{NAD}^+ \\rightarrow \\text{CO}_2 + \\text{NADH}

The reaction model comes from an EnzymeML document
(:data:`FDH_ENZYMEML`), the chemicals are resolved from its identifiers,
and the flowsheet can be exchanged as a DEXPI process model.

"""
from __future__ import annotations

import os
from typing import Optional

__all__ = (
    'FDH_DEXPI',
    'FDH_ENZYMEML',
    'FDH_PROPERTIES',
    'create_fdh_chemicals',
    'create_fdh_system',
)

#: Path of the EnzymeML document of the FDH reaction model.
FDH_ENZYMEML = os.path.join(
    os.path.dirname(__file__), 'data', 'fdh_nadh_regeneration.json'
)

#: Path of a DEXPI process model of the FDH flowsheet as a flowsheeting
#: tool other than BioSTEAM could write it: a feed pump and a reactor,
#: without reactions, residence time or any BioSTEAM annotation.
FDH_DEXPI = os.path.join(
    os.path.dirname(__file__), 'data', 'fdh_flowsheet.xml'
)

#: Properties of the species that chemical databases do not provide: the
#: formula of NADH and the molecular weight of the FDH monomer of
#: *Candida boidinii* (approximate). Remaining properties default to those
#: of water.
FDH_PROPERTIES = {
    'nadh': {'formula': 'C21H29N7O14P2'},
    'fdh': {'MW': 40400.},
}


def create_fdh_chemicals(document=None):
    """
    Return the chemicals of the FDH reference case.

    Parameters
    ----------
    document : EnzymeMLDocument, optional
        Reaction model. Defaults to the document at :data:`FDH_ENZYMEML`.

    Examples
    --------
    >>> from biosteam.interop import examples
    >>> chemicals = examples.create_fdh_chemicals()
    >>> chemicals.IDs
    ('Water', 'formate', 'nad', 'co2', 'nadh', 'fdh')
    >>> chemicals.formate.CAS
    '64-18-6'

    """
    import warnings
    from . import enzymeml
    from .chemicals import PlaceholderWarning
    if document is None:
        document = enzymeml.read(FDH_ENZYMEML)
    with warnings.catch_warnings():
        # NAD+ is found by its InChIKey but lacks property models, which
        # default to those of water; this is intended here.
        warnings.simplefilter('ignore', PlaceholderWarning)
        return document.create_chemicals(properties=FDH_PROPERTIES)


def create_fdh_system(ID: str = 'fdh_sys', *,
                      residence_time: float = 2.,
                      temperature: Optional[float] = None,
                      formate: float = 100.,
                      nad: float = 5.,
                      fdh: float = 2e-4,
                      flow: float = 1000.):
    """
    Create the FDH reference system: a feed pump, a heater and a
    continuous stirred tank reactor in which FDH regenerates NADH.

    The active property package must define the chemicals of
    :func:`create_fdh_chemicals`.

    Parameters
    ----------
    ID :
        ID of the system.
    residence_time :
        Residence time of the reactor [hr].
    temperature :
        Reactor temperature [K]. Defaults to the temperature of the
        EnzymeML document.
    formate, nad, fdh :
        Feed concentrations of formate, NAD+ and enzyme [mol/m3].
    flow :
        Feed flow of water [kg/hr].

    Returns
    -------
    System

    Examples
    --------
    >>> import biosteam as bst
    >>> from biosteam.interop import examples
    >>> bst.main_flowsheet.set_flowsheet('fdh_example')
    >>> bst.settings.set_thermo(examples.create_fdh_chemicals())
    >>> fdh_sys = examples.create_fdh_system()
    >>> fdh_sys.simulate()
    >>> product = fdh_sys.flowsheet.stream.product
    >>> conversion = 1 - (
    ...     product.imol['nad'] / fdh_sys.flowsheet.stream.feed.imol['nad'])
    >>> round(conversion, 3)
    0.912

    """
    import biosteam as bst
    from . import enzymeml
    document = enzymeml.read(FDH_ENZYMEML)
    if temperature is None:
        temperature = document.temperature or 303.15
    volume = flow / 1000.  # Approximate feed volume [m3/hr]
    feed = bst.Stream(
        'feed',
        Water=flow,
        units='kg/hr',
        T=298.15,
    )
    feed.imol['formate'] = formate * volume / 1000.
    feed.imol['nad'] = nad * volume / 1000.
    feed.imol['fdh'] = fdh * volume / 1000.
    P1 = bst.Pump('P1', ins=feed, outs='pressurized_feed', P=2 * 101325)
    H1 = bst.HXutility('H1', ins=P1-0, outs='heated_feed', T=temperature)
    R1 = bst.SinglePhaseReactor(
        'R1', ins=H1-0, outs='product',
        tau=residence_time, T=temperature, P=101325,
    )
    document.parameterize(R1)
    return bst.System.from_units(ID, [P1, H1, R1])
