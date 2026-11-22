# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Reading of EnzymeML documents.

EnzymeML is an exchange format for enzyme kinetics: it records the species
of an enzymatic reaction, its stoichiometry, the rate law and the estimated
kinetic parameters. This module reads the parts of an EnzymeML (version 2)
document that describe the reaction model into immutable records.

"""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field, replace
from typing import Any, Mapping, Optional, Union

from ..chemicals import ChemicalIdentity, is_cas_number, normalize_chebi
from .units import conversion_factor, to_kelvin, unit_name

__all__ = (
    'EnzymeMLDocument',
    'EnzymeMLError',
    'Parameter',
    'Reaction',
    'ReactionElement',
    'Species',
    'read',
)

_CAS_IN_REFERENCE = re.compile(
    r'(?:cas[_:=/ ]*(?:rn=)?|/)(\d{2,7}-\d{2}-\d)\b', re.IGNORECASE)
_PUBCHEM_IN_REFERENCE = re.compile(
    r'pubchem[^\s]*?/compound/(\d+)', re.IGNORECASE)

#: Average masses of amino acid residues [g/mol].
_RESIDUE_MASSES = {
    'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886, 'C': 103.1388,
    'E': 129.1155, 'Q': 128.1307, 'G': 57.0519, 'H': 137.1411, 'I': 113.1594,
    'L': 113.1594, 'K': 128.1741, 'M': 131.1926, 'F': 147.1766, 'P': 97.1167,
    'S': 87.0782, 'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326,
}
_WATER_MASS = 18.01528


class EnzymeMLError(ValueError):
    """Raised when an EnzymeML document cannot be interpreted."""


def protein_molecular_weight(sequence: Optional[str]) -> Optional[float]:
    """
    Return the average molecular weight of a protein [g/mol] from its
    amino acid sequence, or None if the sequence is missing or holds
    non-standard residues.

    Examples
    --------
    >>> from biosteam.interop.enzymeml.document import protein_molecular_weight
    >>> round(protein_molecular_weight('GAG'), 2)
    203.2

    """
    if not sequence:
        return None
    sequence = ''.join(sequence.split()).upper()
    try:
        return sum(_RESIDUE_MASSES[i] for i in sequence) + _WATER_MASS
    except KeyError:
        return None


@dataclass(frozen=True)
class Species:
    """
    A small molecule, protein or complex of an EnzymeML document.

    Attributes
    ----------
    id :
        Identifier used in reactions and rate laws.
    name :
        Name of the species.
    kind :
        Either 'small_molecule', 'protein' or 'complex'.
    inchikey, inchi, smiles :
        Structure identifiers of a small molecule.
    references :
        Links to database entries (e.g. ChEBI, PubChem, UniProt).
    ecnumber, organism, sequence :
        Enzyme classification, source organism and amino acid sequence of
        a protein.

    """

    id: str
    name: str = ''
    kind: str = 'small_molecule'
    inchikey: Optional[str] = None
    inchi: Optional[str] = None
    smiles: Optional[str] = None
    references: tuple[str, ...] = ()
    ecnumber: Optional[str] = None
    organism: Optional[str] = None
    sequence: Optional[str] = None

    def _reference_match(self, pattern) -> Optional[str]:
        for reference in self.references:
            match = pattern.search(reference)
            if match:
                return match.group(1)
        return None

    def to_identity(self, ID: Optional[str] = None) -> ChemicalIdentity:
        """
        Return the chemical identity of the species.

        Parameters
        ----------
        ID :
            ID of the chemical in BioSTEAM. Defaults to the species id.

        """
        chebi = None
        for reference in self.references:
            chebi = (
                normalize_chebi(reference) if 'chebi' in reference.lower()
                else None
            )
            if chebi:
                break
        CAS = self._reference_match(_CAS_IN_REFERENCE)
        is_small_molecule = self.kind == 'small_molecule'
        return ChemicalIdentity(
            ID=ID or self.id,
            name=self.name or None,
            CAS=CAS if is_cas_number(CAS) else None,
            InChIKey=self.inchikey,
            InChI=self.inchi,
            SMILES=self.smiles,
            ChEBI=chebi,
            PubChem=self._reference_match(_PUBCHEM_IN_REFERENCE),
            MW=(None if is_small_molecule
                else protein_molecular_weight(self.sequence)),
            phase=None if is_small_molecule else 'l',
            pure=is_small_molecule,
        )


@dataclass(frozen=True)
class Parameter:
    """
    A kinetic parameter of an EnzymeML document.

    Attributes
    ----------
    id, symbol, name :
        Identifier, symbol used in rate laws, and name.
    value :
        Value in the unit of the document; None if not estimated.
    unit :
        EnzymeML unit definition of the value.

    """

    id: str
    symbol: str
    name: str = ''
    value: Optional[float] = None
    unit: Optional[Mapping[str, Any]] = field(default=None, compare=False)

    @property
    def unit_name(self) -> Optional[str]:
        """Display name of the unit."""
        return unit_name(self.unit)

    @property
    def coherent_value(self) -> Optional[float]:
        """Value in the coherent system of units (kmol, m, hr, kg, K)."""
        if self.value is None:
            return None
        return self.value * conversion_factor(self.unit)


@dataclass(frozen=True)
class ReactionElement:
    """A reactant or product with its (positive) stoichiometric number."""

    species_id: str
    stoichiometry: float = 1.


@dataclass(frozen=True)
class Reaction:
    """
    A reaction of an EnzymeML document.

    Attributes
    ----------
    id, name :
        Identifier and name.
    reversible :
        Whether the reaction is declared reversible.
    reactants, products :
        Species consumed and formed.
    modifiers :
        Role of each species that takes part without being converted, e.g.
        ``{'fdh': 'biocatalyst'}``.
    rate_law :
        Rate law expression in terms of species ids and parameter symbols.

    """

    id: str
    name: str = ''
    reversible: bool = False
    reactants: tuple[ReactionElement, ...] = ()
    products: tuple[ReactionElement, ...] = ()
    modifiers: Mapping[str, str] = field(default_factory=dict)
    rate_law: Optional[str] = None

    @property
    def stoichiometry(self) -> dict[str, float]:
        """Signed stoichiometric number of each species."""
        stoichiometry: dict[str, float] = {}
        for element in self.reactants:
            stoichiometry[element.species_id] = (
                stoichiometry.get(element.species_id, 0.)
                - element.stoichiometry
            )
        for element in self.products:
            stoichiometry[element.species_id] = (
                stoichiometry.get(element.species_id, 0.)
                + element.stoichiometry
            )
        return stoichiometry


@dataclass(frozen=True)
class EnzymeMLDocument:
    """
    The reaction model of an EnzymeML document.

    Create documents with :func:`read`.

    Attributes
    ----------
    name :
        Title of the document.
    species :
        Small molecules, proteins and complexes by id.
    reactions :
        Reactions by id.
    parameters :
        Kinetic parameters by symbol.
    temperature :
        Temperature of the measurements [K], if they agree on one.
    pH :
        pH of the measurements, if they agree on one.
    chemical_IDs :
        ID of the BioSTEAM chemical for each species id. Species that are
        not listed use their id.

    Examples
    --------
    Load the reference data set of the NADH regeneration by formate
    dehydrogenase and create the chemicals and the kinetic reaction:

    >>> import biosteam as bst
    >>> from biosteam.interop import enzymeml, examples
    >>> document = enzymeml.read(examples.FDH_ENZYMEML)
    >>> sorted(document.species)
    ['co2', 'fdh', 'formate', 'nad', 'nadh']
    >>> document.reactions['r_fdh'].rate_law
    'kcat * fdh * nad * formate / (K_i_nad * K_m_formate + K_m_formate * nad
    + K_m_nad * formate + nad * formate)'
    >>> (document.parameters['kcat'].value,
    ...  document.parameters['kcat'].unit_name)
    (3.7, '1 / s')

    """

    name: str = ''
    species: Mapping[str, Species] = field(default_factory=dict)
    reactions: Mapping[str, Reaction] = field(default_factory=dict)
    parameters: Mapping[str, Parameter] = field(default_factory=dict)
    temperature: Optional[float] = None
    pH: Optional[float] = None
    chemical_IDs: Mapping[str, str] = field(default_factory=dict)

    def with_chemical_IDs(
            self, chemical_IDs: Mapping[str, str]) -> EnzymeMLDocument:
        """
        Return a copy of the document that maps species to other chemical
        IDs, e.g. ``document.with_chemical_IDs({'s0': 'Formate'})``.
        """
        unknown = sorted(set(chemical_IDs).difference(self.species))
        if unknown:
            raise EnzymeMLError(f"document has no species {unknown}")
        return replace(
            self, chemical_IDs={**self.chemical_IDs, **chemical_IDs})

    def chemical_ID(self, species_id: str) -> str:
        """Return the ID of the BioSTEAM chemical for a species."""
        return self.chemical_IDs.get(species_id, species_id)

    def get_reaction(self, reaction_id: Optional[str] = None) -> Reaction:
        """
        Return a reaction by id; the id may be omitted if the document
        has a single reaction.
        """
        if reaction_id is None:
            if len(self.reactions) != 1:
                raise EnzymeMLError(
                    f"document has {len(self.reactions)} reactions; "
                    f"choose one of {sorted(self.reactions)}"
                )
            return next(iter(self.reactions.values()))
        try:
            return self.reactions[reaction_id]
        except KeyError:
            raise EnzymeMLError(
                f"document has no reaction '{reaction_id}'; "
                f"choose one of {sorted(self.reactions)}"
            ) from None

    def chemical_identities(self) -> list[ChemicalIdentity]:
        """Return the chemical identity of every species."""
        return [
            species.to_identity(self.chemical_ID(species.id))
            for species in self.species.values()
        ]

    def create_chemicals(self, resolver=None, *, chemicals=None,
                         properties=None, online=False, extra=()):
        """
        Return the chemicals of all species, with water added.

        Parameters
        ----------
        resolver : ResolverChain, optional
            Resolver of chemical identities. Defaults to
            :func:`~biosteam.interop.chemicals.default_resolver` built
            from the remaining arguments.
        chemicals : Iterable[Chemical], optional
            Chemicals defined by the user; these take precedence.
        properties : Mapping, optional
            Property records for species that databases do not know, such
            as enzymes, e.g. ``{'fdh': {'MW': 40370.}}``.
        online : bool, optional
            Whether to query PubChem for unknown species.
        extra : Iterable[str], optional
            IDs of further chemicals to include. Defaults to no further
            chemicals; water is always included as the solvent.

        Returns
        -------
        Chemicals
            Compiled chemicals, ready for ``settings.set_thermo``.

        """
        from .kinetics import create_chemicals
        return create_chemicals(
            self, resolver, chemicals=chemicals, properties=properties,
            online=online, extra=extra,
        )

    def to_reaction(self, reaction_id: Optional[str] = None, *,
                    residence_time: Optional[float] = None,
                    reactant: Optional[str] = None,
                    chemicals=None,
                    flow_pattern: str = 'mixed'):
        """
        Return a reaction of the document as a kinetic
        :class:`thermosteam.Reaction` for continuous stirred tank reactors.

        Parameters
        ----------
        reaction_id :
            Id of the reaction. May be omitted if there is only one.
        residence_time :
            Residence time of the reactor [hr]. Alternatively, bind the
            reaction to a reactor with :meth:`parameterize`.
        reactant :
            Species id of the reactant that conversion refers to. Defaults
            to the first reactant.
        chemicals : Chemicals, optional
            Chemicals of the reaction. Defaults to the active property
            package.
        flow_pattern :
            Flow pattern of the reactor: 'mixed' for an ideal stirred tank
            or 'plug' for an ideal tubular reactor.

        """
        from .kinetics import create_reaction
        return create_reaction(
            self, reaction_id, residence_time=residence_time,
            reactant=reactant, chemicals=chemicals,
            flow_pattern=flow_pattern,
        )

    def parameterize(self, reactor, reaction_id: Optional[str] = None, *,
                     reactant: Optional[str] = None,
                     flow_pattern: Optional[str] = None):
        """
        Set the reaction of a stirred tank reactor from this document.

        The kinetics follow the residence time of the reactor (its `tau`
        attribute), also when it is changed later.

        Parameters
        ----------
        reactor : Unit
            A reactor with `reactions` and `tau` attributes, such as
            :class:`~biosteam.units.SinglePhaseReactor`.
        reaction_id :
            Id of the reaction. May be omitted if there is only one.
        reactant :
            Species id of the reactant that conversion refers to.
        flow_pattern :
            'mixed' (ideal stirred tank) or 'plug' (ideal tubular
            reactor). Defaults to the `flow_pattern` attribute of the
            reactor, which is set when a reactor is imported from DEXPI,
            and otherwise to 'mixed'.

        Returns
        -------
        Reaction
            The reaction assigned to the reactor.

        """
        if flow_pattern is None:
            flow_pattern = getattr(reactor, 'flow_pattern', None) or 'mixed'
        reaction = self.to_reaction(
            reaction_id, reactant=reactant, chemicals=reactor.chemicals,
            flow_pattern=flow_pattern,
        )
        reaction.rate.bind(reactor)
        reactor.reactions = reaction
        return reaction

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> EnzymeMLDocument:
        """Create a document from the content of an EnzymeML JSON file."""
        if not isinstance(data, Mapping):
            raise EnzymeMLError("an EnzymeML document must be a JSON object")
        species: dict[str, Species] = {}
        for kind, key in (('small_molecule', 'small_molecules'),
                          ('protein', 'proteins'), ('complex', 'complexes')):
            for item in data.get(key) or ():
                one = _species(item, kind)
                if one.id in species:
                    raise EnzymeMLError(f"duplicate species id '{one.id}'")
                species[one.id] = one
        parameters: dict[str, Parameter] = {}
        for item in data.get('parameters') or ():
            parameter = Parameter(
                id=_required(item, 'id', 'parameter'),
                symbol=item.get('symbol') or item['id'],
                name=item.get('name') or '',
                value=_number(item.get('value')),
                unit=item.get('unit'),
            )
            parameters[parameter.symbol] = parameter
        rate_laws = {
            item.get('species_id'): item.get('equation')
            for item in data.get('equations') or ()
            if _text(item.get('equation_type')).lower() == 'ratelaw'
        }
        reactions: dict[str, Reaction] = {}
        for item in data.get('reactions') or ():
            reaction = _reaction(item, rate_laws)
            for element in (*reaction.reactants, *reaction.products):
                if element.species_id not in species:
                    raise EnzymeMLError(
                        f"reaction '{reaction.id}' refers to undefined "
                        f"species '{element.species_id}'"
                    )
            reactions[reaction.id] = reaction
        measurements = data.get('measurements') or ()
        temperatures = {
            to_kelvin(_number(i.get('temperature')), i.get('temperature_unit'))
            for i in measurements
        }
        pHs = {_number(i.get('ph')) for i in measurements}
        return cls(
            name=data.get('name') or '',
            species=species,
            reactions=reactions,
            parameters=parameters,
            temperature=temperatures.pop() if len(temperatures) == 1 else None,
            pH=pHs.pop() if len(pHs) == 1 else None,
        )


def _text(value: Any) -> str:
    return str(getattr(value, 'value', value) or '')


def _number(value: Any) -> Optional[float]:
    if value is None or isinstance(value, bool):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        raise EnzymeMLError(f"expected a number, found {value!r}") from None


def _required(item: Mapping, key: str, what: str) -> str:
    value = item.get(key) if isinstance(item, Mapping) else None
    if not value:
        raise EnzymeMLError(f"{what} without '{key}'")
    return str(value)


def _species(item: Mapping, kind: str) -> Species:
    return Species(
        id=_required(item, 'id', kind),
        name=item.get('name') or '',
        kind=kind,
        inchikey=item.get('inchikey') or None,
        inchi=item.get('inchi') or None,
        smiles=item.get('canonical_smiles') or None,
        references=tuple(item.get('references') or ()),
        ecnumber=item.get('ecnumber') or None,
        organism=item.get('organism') or None,
        sequence=item.get('sequence') or None,
    )


def _elements(items) -> tuple[ReactionElement, ...]:
    elements = []
    for item in items or ():
        stoichiometry = _number(item.get('stoichiometry'))
        elements.append(ReactionElement(
            _required(item, 'species_id', 'reaction element'),
            abs(stoichiometry) if stoichiometry else 1.,
        ))
    return tuple(elements)


def _reaction(item: Mapping, rate_laws: Mapping) -> Reaction:
    reaction_id = _required(item, 'id', 'reaction')
    kinetic_law = item.get('kinetic_law') or {}
    return Reaction(
        id=reaction_id,
        name=item.get('name') or '',
        reversible=bool(item.get('reversible')),
        reactants=_elements(item.get('reactants')),
        products=_elements(item.get('products')),
        modifiers={
            _required(i, 'species_id', 'modifier'): _text(i.get('role'))
            for i in item.get('modifiers') or ()
        },
        rate_law=kinetic_law.get('equation') or rate_laws.get(reaction_id),
    )


def _read_omex(path: str) -> Mapping[str, Any]:
    """Read an OMEX archive with pyenzyme."""
    try:
        import pyenzyme
    except ImportError:
        raise ImportError(
            "reading EnzymeML OMEX archives requires the 'pyenzyme' "
            "package; install it with `pip install biosteam[enzymeml]` or "
            "export the document as EnzymeML JSON"
        ) from None
    return pyenzyme.read_enzymeml(path).model_dump(mode='json')


def read(
        source: Union[str, os.PathLike, Mapping[str, Any], Any]
) -> EnzymeMLDocument:
    """
    Read an EnzymeML document.

    Parameters
    ----------
    source :
        One of the following:

        * path of an EnzymeML JSON file (``.json``);
        * path of an EnzymeML OMEX archive (``.omex``), which requires the
          optional `pyenzyme` package;
        * the content of an EnzymeML JSON file as a dictionary;
        * a ``pyenzyme.EnzymeMLDocument`` object.

    Returns
    -------
    EnzymeMLDocument

    Raises
    ------
    EnzymeMLError
        If the document is not valid EnzymeML.
    ImportError
        If an OMEX archive is given and `pyenzyme` is not installed.

    """
    if isinstance(source, EnzymeMLDocument):
        return source
    if isinstance(source, Mapping):
        return EnzymeMLDocument.from_dict(source)
    if hasattr(source, 'model_dump'):  # pyenzyme (pydantic) document
        return EnzymeMLDocument.from_dict(source.model_dump(mode='json'))
    path = os.fspath(source)
    if path.lower().endswith('.omex'):
        return EnzymeMLDocument.from_dict(_read_omex(path))
    with open(path, encoding='utf-8') as file:
        try:
            data = json.load(file)
        except ValueError as error:
            raise EnzymeMLError(
                f"'{path}' is not an EnzymeML JSON document: {error}"
            ) from None
    return EnzymeMLDocument.from_dict(data)
