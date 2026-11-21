# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
A stand-in unit operation for process steps without a BioSTEAM model.

"""
from __future__ import annotations

from typing import Optional, Sequence

import thermosteam as tmo

from ..._unit import Unit

__all__ = ('GenericProcessStep',)


class GenericProcessStep(Unit):
    """
    Create a placeholder for a DEXPI process step that has no equivalent
    BioSTEAM unit operation.

    The unit keeps the flowsheet connected and solvable: it mixes all
    inlets and distributes the mixture over the outlets. It performs no
    design or costing. The original DEXPI class and attributes are kept and
    written back when the system is exported again.

    Parameters
    ----------
    ins :
        Inlet streams.
    outs :
        Outlet streams.
    dexpi_class :
        Name of the DEXPI process step class this unit stands for.
    split :
        Fraction of the mixed inlets sent to each outlet. Defaults to an
        even distribution.

    Examples
    --------
    >>> import biosteam as bst
    >>> from biosteam.interop.dexpi import GenericProcessStep
    >>> bst.main_flowsheet.set_flowsheet('generic_step_example')
    >>> bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> feed = bst.Stream('feed', Water=10, Ethanol=2)
    >>> step = GenericProcessStep(
    ...     'U1', ins=feed, outs=('a', 'b'), dexpi_class='Sieving',
    ...     split=(0.25, 0.75))
    >>> step.simulate()
    >>> step.outs[0].imol['Water'], step.outs[1].imol['Water']
    (2.5, 7.5)

    """

    _N_ins = 1
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    line = 'Generic process step'

    def _init(self, dexpi_class: str = 'ProcessStep',
              split: Optional[Sequence[float]] = None):
        #: Name of the DEXPI class this unit stands for.
        self.dexpi_class = dexpi_class
        #: Fraction of the mixed inlets sent to each outlet.
        self.split = None if split is None else tuple(split)
        #: The imported DEXPI object, if any; its attributes are written
        #: back on export.
        self.dexpi_step = None
        self._mixture = tmo.Stream(None, thermo=self.thermo)

    def _run(self):
        outs = [i for i in self.outs if i]
        ins = [i for i in self.ins if i]
        if not outs or not ins:
            return
        mixture = self._mixture
        mixture.mix_from(ins)
        split = self.split
        if split is None:
            split = (1. / len(outs),) * len(outs)
        elif len(split) != len(outs):
            raise ValueError(
                f"length of split ({len(split)}) must match the number of "
                f"outlets ({len(outs)})"
            )
        mol = mixture.mol
        for outlet, fraction in zip(outs, split):
            outlet.copy_like(mixture)
            outlet.mol[:] = mol * fraction
