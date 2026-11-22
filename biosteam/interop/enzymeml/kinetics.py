# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Kinetic reactions for continuous stirred tank reactors from EnzymeML.

"""
from __future__ import annotations

from typing import Any, Mapping, Optional

import thermosteam as tmo
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from ..chemicals import ChemicalIdentity, default_resolver
from .document import EnzymeMLDocument, EnzymeMLError
from .ratelaw import RateLaw

__all__ = (
    'PlugFlowRate',
    'SteadyStateCSTRRate',
    'create_chemicals',
    'create_reaction',
)

_WATER = ChemicalIdentity('Water', CAS='7732-18-5', name='water')


class SteadyStateCSTRRate:
    """
    The rate of a kinetically controlled reaction in a continuous stirred
    tank reactor (CSTR) at steady state.

    A :class:`thermosteam.Reaction` with a `rate` asks for the amount of
    reactant converted per hour, given the mixed feed. In an ideal CSTR the
    reaction proceeds at outlet conditions, so the conversion :math:`\\xi`
    [kmol/hr] is the solution of the design equation

    .. math:: \\xi = r(C_{out}(\\xi)) \\cdot V

    with outlet concentrations :math:`C_{out,i} = (F_i + \\nu_i \\xi) / Q`
    for molar feed flows :math:`F_i`, stoichiometric numbers :math:`\\nu_i`
    (relative to the reactant) and volumetric flow :math:`Q`, and with the
    liquid volume :math:`V = \\tau Q` for residence time :math:`\\tau`. The
    volumetric flow is assumed not to change by reaction.

    Parameters
    ----------
    rate_law :
        Rate of reaction [kmol/m3/hr] as an expression of species
        concentrations [kmol/m3] and parameters.
    parameters :
        Value of each parameter of the rate law, in units coherent with
        kmol, m3 and hr.
    species :
        ID of the chemical for each species name used in the rate law.
    stoichiometry :
        Stoichiometric number of each chemical (negative for reactants).
    reactant :
        ID of the reactant that the conversion refers to.
    residence_time :
        Residence time of the reactor [hr]. Not needed if the rate is bound
        to a reactor with :meth:`bind`.
    source :
        Provenance of the kinetics (free-form, JSON-compatible).

    Examples
    --------
    A first-order reaction A -> B with k = 2 1/hr in a CSTR with a
    residence time of 0.5 hr converts k*tau / (1 + k*tau) = 50% of A:

    >>> import biosteam as bst
    >>> from biosteam.interop.enzymeml import SteadyStateCSTRRate
    >>> bst.main_flowsheet.set_flowsheet('cstr_rate_example')
    >>> bst.settings.set_thermo(['Water', 'Ethanol', 'Propanol'], cache=True)
    >>> rate = SteadyStateCSTRRate(
    ...     'k * A', {'k': 2.}, species={'A': 'Ethanol'},
    ...     stoichiometry={'Ethanol': -1, 'Propanol': 1}, reactant='Ethanol')
    >>> reaction = bst.Reaction('Ethanol -> Propanol', rate=rate)
    >>> feed = bst.Stream('feed', Water=100, Ethanol=1)
    >>> R1 = bst.SinglePhaseReactor('R1', ins=feed, reactions=reaction,
    ...                             tau=0.5, T=298.15, P=101325)
    >>> rate.bind(R1)
    >>> R1.run()
    >>> round(R1.outs[0].imol['Ethanol'], 6)
    0.5

    """

    type_name = 'SteadyStateCSTRRate'

    def __init__(self, rate_law, parameters: Mapping[str, float],
                 species: Mapping[str, str],
                 stoichiometry: Mapping[str, float],
                 reactant: str,
                 residence_time: Optional[float] = None,
                 source: Optional[Mapping[str, Any]] = None):
        self.parameters = {i: float(j) for i, j in parameters.items()}
        self.species = dict(species)
        names = set(self.parameters) | set(self.species)
        if not isinstance(rate_law, RateLaw):
            rate_law = RateLaw(rate_law, names)
        self.rate_law = rate_law
        self.stoichiometry = {
            i: float(j) for i, j in stoichiometry.items() if j
        }
        if self.stoichiometry.get(reactant, 0.) >= 0.:
            raise ValueError(f"'{reactant}' is not a reactant of the reaction")
        self.reactant = reactant
        self.source = dict(source) if source else None
        self._residence_time = residence_time
        self._reactor = None
        # Conversion is expressed per unit of reactant consumed, whereas
        # the rate law counts reaction events.
        self._reactant_number = -self.stoichiometry[reactant]

    # %% Reactor volume

    def bind(self, reactor) -> None:
        """
        Use the residence time of `reactor` (its `tau` attribute), also
        when it changes after binding.
        """
        if not hasattr(reactor, 'tau'):
            raise TypeError(
                f"'{type(reactor).__name__}' object has no residence time "
                "attribute 'tau'"
            )
        self._reactor = reactor

    @property
    def residence_time(self) -> Optional[float]:
        """Residence time of the reactor [hr]."""
        if self._reactor is not None:
            return self._reactor.tau
        return self._residence_time

    @residence_time.setter
    def residence_time(self, residence_time):
        self._reactor = None
        self._residence_time = residence_time

    # %% Rate

    def volumetric_rate(self, concentrations: Mapping[str, float]) -> float:
        """
        Return the rate of reaction [kmol/m3/hr] per unit of reactant.

        Parameters
        ----------
        concentrations :
            Concentration [kmol/m3] of each chemical in the rate law by ID.

        """
        values = dict(self.parameters)
        for name, ID in self.species.items():
            values[name] = concentrations[ID]
        return self.rate_law(values) * self._reactant_number

    def __call__(self, stream) -> float:
        """
        Return the conversion of the reactant [kmol/hr] at steady state.

        Parameters
        ----------
        stream : Stream
            The mixed feed of the reactor at operating conditions.

        """
        tau = self.residence_time
        if tau is None:
            raise RuntimeError(
                "residence time is not defined; set `residence_time` or "
                "bind the rate to a reactor"
            )
        flow = stream.F_vol  # m3/hr
        if flow <= 0. or tau <= 0.:
            return 0.
        IDs = set(self.species.values()) | set(self.stoichiometry)
        imol = stream.imol
        feed = {ID: float(imol[ID]) for ID in IDs}
        reactant_number = self._reactant_number
        changes = {
            ID: nu / reactant_number for ID, nu in self.stoichiometry.items()
        }

        def rate(extent):
            """Rate [kmol/m3/hr] after converting `extent` [kmol/hr]."""
            return self.volumetric_rate({
                ID: max(feed[ID] + changes.get(ID, 0.) * extent, 0.) / flow
                for ID in IDs
            })

        forward_limit = min(
            feed[ID] / -change for ID, change in changes.items() if change < 0.
        )
        reverse_limit = -min(
            (feed[ID] / change for ID, change in changes.items()
             if change > 0.),
            default=0.,
        )
        return self._solve(rate, tau, flow, reverse_limit, forward_limit)

    def _solve(self, rate, tau, flow, lower_limit, upper_limit) -> float:
        """Return the conversion [kmol/hr] for a rate as a function of
        conversion, within the limits set by the feed."""
        volume = tau * flow

        def residual(extent):
            return rate(extent) * volume - extent

        at_feed = residual(0.)
        if at_feed > 0.:
            lower, upper = 0., upper_limit
            if upper <= 0.:
                return 0.
            if residual(upper) >= 0.:
                return upper  # Reaction is limited by the reactants
        elif at_feed < 0.:
            lower, upper = lower_limit, 0.
            if lower >= 0.:
                return 0.
            if residual(lower) <= 0.:
                return lower  # Reverse reaction is limited by the products
        else:
            return 0.
        return brentq(
            residual, lower, upper,
            xtol=1e-15 * (upper - lower), rtol=1e-15, maxiter=200,
        )

    # %% Serialization

    def to_dict(self) -> dict:
        """Return the rate as a dictionary of JSON-compatible values."""
        data = {
            'type': self.type_name,
            'rate_law': self.rate_law.equation,
            'units': 'kmol, m3, hr',
            'parameters': dict(self.parameters),
            'species': dict(self.species),
            'stoichiometry': dict(self.stoichiometry),
            'reactant': self.reactant,
        }
        if self._reactor is None and self._residence_time is not None:
            data['residence_time'] = float(self._residence_time)
        if self.source:
            data['source'] = dict(self.source)
        return data

    @classmethod
    def from_dict(cls, data: Mapping[str, Any],
                  chemicals=None) -> SteadyStateCSTRRate:
        """Create a rate from a dictionary written by :meth:`to_dict`."""
        return cls(
            data['rate_law'],
            data['parameters'],
            data['species'],
            data['stoichiometry'],
            data['reactant'],
            residence_time=data.get('residence_time'),
            source=data.get('source'),
        )

    def __repr__(self):
        return f"{type(self).__name__}({self.rate_law.equation!r})"


class PlugFlowRate(SteadyStateCSTRRate):
    """
    The rate of a kinetically controlled reaction in an ideal plug-flow
    (tubular) reactor at steady state.

    The conversion :math:`\\xi` [kmol/hr] follows from integrating the
    rate along the reactor,

    .. math:: \\frac{d\\xi}{dt} = r(C(\\xi)) \\cdot Q, \\qquad \\xi(0) = 0,

    up to the residence time :math:`\\tau`, with concentrations
    :math:`C_i = (F_i + \\nu_i \\xi) / Q`. The volumetric flow :math:`Q`
    is assumed not to change by reaction. Parameters are those of
    :class:`SteadyStateCSTRRate`.

    Examples
    --------
    A first-order reaction A -> B with k = 2 1/hr converts
    1 - exp(-k*tau) = 63.2% of A at a residence time of 0.5 hr:

    >>> import biosteam as bst
    >>> from biosteam.interop.enzymeml import PlugFlowRate
    >>> bst.main_flowsheet.set_flowsheet('plug_flow_rate_example')
    >>> bst.settings.set_thermo(['Water', 'Ethanol', 'Propanol'], cache=True)
    >>> rate = PlugFlowRate(
    ...     'k * A', {'k': 2.}, species={'A': 'Ethanol'},
    ...     stoichiometry={'Ethanol': -1, 'Propanol': 1},
    ...     reactant='Ethanol', residence_time=0.5)
    >>> feed = bst.Stream('feed', Water=100, Ethanol=1)
    >>> round(rate(feed), 6)
    0.632121

    """

    type_name = 'PlugFlowRate'

    def _solve(self, rate, tau, flow, lower_limit, upper_limit) -> float:
        at_feed = rate(0.)
        limit = upper_limit if at_feed > 0. else lower_limit
        if at_feed == 0. or limit == 0.:
            return 0.

        # Integrate the fraction of the limiting conversion, which is of
        # order one whatever the scale of the flows.
        def derivative(time, fraction):
            return [rate(fraction[0] * limit) * flow / limit]

        solution = solve_ivp(
            derivative, (0., tau), [0.], method='LSODA',
            rtol=1e-10, atol=1e-13,
        )
        if not solution.success:
            raise RuntimeError(
                f"integration of the rate law failed: {solution.message}"
            )
        return min(max(solution.y[0, -1], 0.), 1.) * limit


#: Rate classes by flow pattern of the reactor.
RATE_CLASSES = {
    'mixed': SteadyStateCSTRRate,
    'plug': PlugFlowRate,
}


def create_chemicals(document: EnzymeMLDocument, resolver=None, *,
                     chemicals=None, properties=None, online=False,
                     extra=()) -> tmo.Chemicals:
    """
    Return the chemicals of all species of a document, with water added.

    See :meth:`EnzymeMLDocument.create_chemicals`.
    """
    if resolver is None:
        resolver = default_resolver(
            chemicals=chemicals, properties=properties, online=online,
        )
    identities, IDs = [], set()
    for identity in document.chemical_identities():
        # Species that share a chemical ID (e.g. the same compound in
        # several vessels) are one chemical.
        if identity.ID not in IDs:
            identities.append(identity)
            IDs.add(identity.ID)
    if 'Water' not in IDs and not any(i.CAS == _WATER.CAS for i in identities):
        identities.insert(0, _WATER)
        IDs.add('Water')
    for ID in extra:
        if ID not in IDs:
            identities.append(ChemicalIdentity(ID))
            IDs.add(ID)
    return resolver.create_chemicals(identities)


def create_reaction(document: EnzymeMLDocument,
                    reaction_id: Optional[str] = None, *,
                    residence_time: Optional[float] = None,
                    reactant: Optional[str] = None,
                    chemicals: Optional[tmo.Chemicals] = None,
                    flow_pattern: str = 'mixed') -> tmo.Reaction:
    """
    Return a reaction of a document as a kinetic reaction for CSTRs.

    See :meth:`EnzymeMLDocument.to_reaction`.
    """
    try:
        rate_class = RATE_CLASSES[flow_pattern]
    except KeyError:
        raise ValueError(
            f"flow_pattern must be one of {sorted(RATE_CLASSES)}, "
            f"not {flow_pattern!r}"
        ) from None
    reaction = document.get_reaction(reaction_id)
    if not reaction.reactants:
        raise EnzymeMLError(f"reaction '{reaction.id}' has no reactants")
    if not reaction.rate_law:
        raise EnzymeMLError(
            f"reaction '{reaction.id}' has no rate law; kinetic parameters "
            "must be estimated before the reaction can be simulated"
        )
    chemical_ID = document.chemical_ID
    if reactant is None:
        reactant = reaction.reactants[0].species_id
    elif reactant not in {i.species_id for i in reaction.reactants}:
        raise EnzymeMLError(
            f"'{reactant}' is not a reactant of reaction '{reaction.id}'"
        )
    names = set(document.species) | set(document.parameters)
    rate_law = RateLaw(reaction.rate_law, names)
    parameters = {}
    for symbol in sorted(rate_law.symbols.intersection(document.parameters)):
        value = document.parameters[symbol].coherent_value
        if value is None:
            raise EnzymeMLError(
                f"parameter '{symbol}' of reaction '{reaction.id}' "
                "has no value"
            )
        parameters[symbol] = value
    species = {
        i: chemical_ID(i)
        for i in sorted(rate_law.symbols.intersection(document.species))
        if i not in parameters
    }
    stoichiometry = {
        chemical_ID(i): j for i, j in reaction.stoichiometry.items()
    }
    enzymes = [
        document.species[i] for i in reaction.modifiers
        if i in document.species
        and document.species[i].kind != 'small_molecule'
    ]
    source = {
        'format': 'EnzymeML',
        'document': document.name,
        'reaction': reaction.id,
    }
    if enzymes:
        source['enzymes'] = [
            {i: j for i, j in (('id', e.id), ('name', e.name),
                               ('ecnumber', e.ecnumber),
                               ('organism', e.organism)) if j}
            for e in enzymes
        ]
    if document.temperature is not None:
        source['temperature'] = document.temperature
    if document.pH is not None:
        source['pH'] = document.pH
    rate = rate_class(
        rate_law, parameters, species, stoichiometry,
        reactant=chemical_ID(reactant),
        residence_time=residence_time, source=source,
    )
    return tmo.Reaction(
        stoichiometry, chemical_ID(reactant), chemicals=chemicals, rate=rate,
    )
