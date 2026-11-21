# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Chemical identities and their resolution to thermosteam chemicals.

Exchange formats such as DEXPI and EnzymeML say *which* chemicals take part
in a process, not how their thermophysical properties are estimated. This
module keeps the two concerns apart:

* :class:`ChemicalIdentity` is a plain record of the identifiers a document
  provides for one chemical;
* a *resolver* turns an identity into a :class:`thermosteam.Chemical`;
* :class:`ResolverChain` tries resolvers in order, so that user data takes
  precedence over databases, and databases over placeholders.

"""
from __future__ import annotations

import json
import re
import urllib.error
import urllib.parse
import urllib.request
import warnings
from dataclasses import asdict, dataclass, field, replace
from typing import Any, Callable, Iterable, Mapping, Optional, Protocol

import thermosteam as tmo

__all__ = (
    'ChemicalIdentity',
    'ChemicalResolutionError',
    'ChemicalResolver',
    'PlaceholderResolver',
    'PlaceholderWarning',
    'PropertyDatabaseResolver',
    'PubChemResolver',
    'ResolverChain',
    'ThermosteamResolver',
    'UserChemicalsResolver',
    'default_resolver',
)

_CAS_PATTERN = re.compile(r'^\d{2,7}-\d{2}-\d$')
_CHEBI_PATTERN = re.compile(r'CHEBI[:_](\d+)', re.IGNORECASE)
_FORMULA_FROM_INCHI = re.compile(r'^InChI=1S?/([^/]+)')

#: Largest molar heat capacity accepted by the property models [J/mol/K].
_MAX_MOLAR_HEAT_CAPACITY = 1e4

#: Heat capacity of water [J/g/K], the default for placeholder chemicals.
_WATER_HEAT_CAPACITY = 4.184


class ChemicalResolutionError(LookupError):
    """Raised when no resolver can provide a chemical."""


class PlaceholderWarning(UserWarning):
    """Warns that a chemical uses default (water-like) properties."""


def is_cas_number(text: Optional[str]) -> bool:
    """
    Return whether `text` is formatted as a CAS registry number.

    Examples
    --------
    >>> from biosteam.interop.chemicals import is_cas_number
    >>> is_cas_number('7732-18-5'), is_cas_number('Water')
    (True, False)

    """
    return bool(text) and bool(_CAS_PATTERN.match(text))


def normalize_chebi(text: Optional[str]) -> Optional[str]:
    """
    Return a ChEBI identifier as ``'CHEBI:<number>'``, or None.

    Examples
    --------
    >>> from biosteam.interop.chemicals import normalize_chebi
    >>> normalize_chebi('https://purl.obolibrary.org/obo/CHEBI_15740')
    'CHEBI:15740'
    >>> normalize_chebi('15740')
    'CHEBI:15740'

    """
    if not text:
        return None
    match = _CHEBI_PATTERN.search(text)
    if match:
        return f'CHEBI:{match.group(1)}'
    if text.strip().isdigit():
        return f'CHEBI:{text.strip()}'
    return None


@dataclass(frozen=True)
class ChemicalIdentity:
    """
    The identifiers a document provides for one chemical.

    Parameters
    ----------
    ID :
        Identifier of the chemical within BioSTEAM.
    name :
        Common name.
    CAS :
        CAS registry number.
    InChIKey :
        IUPAC international chemical identifier key.
    InChI :
        IUPAC international chemical identifier.
    SMILES :
        Simplified molecular-input line-entry specification.
    ChEBI :
        ChEBI identifier, as ``'CHEBI:<number>'``.
    PubChem :
        PubChem compound identifier.
    IUPAC_name :
        Systematic name.
    formula :
        Chemical formula.
    MW :
        Molecular weight [g/mol].
    phase :
        Phase of the chemical if it does not take part in phase equilibrium
        ('s', 'l' or 'g').
    pure :
        Whether the chemical is a defined compound, as opposed to a
        pseudo-component such as an enzyme preparation.

    Examples
    --------
    >>> from biosteam.interop.chemicals import ChemicalIdentity
    >>> water = ChemicalIdentity('H2O', CAS='7732-18-5', name='water')
    >>> water.search_keys()
    ('7732-18-5', 'water', 'H2O')

    """

    ID: str
    name: Optional[str] = None
    CAS: Optional[str] = None
    InChIKey: Optional[str] = None
    InChI: Optional[str] = None
    SMILES: Optional[str] = None
    ChEBI: Optional[str] = None
    PubChem: Optional[str] = None
    IUPAC_name: Optional[str] = None
    formula: Optional[str] = None
    MW: Optional[float] = None
    phase: Optional[str] = None
    pure: bool = True

    def search_keys(self) -> tuple[str, ...]:
        """Return database search terms, most specific first."""
        keys = []
        if is_cas_number(self.CAS):
            keys.append(self.CAS)
        if self.InChIKey:
            keys.append(f'InChIKey={self.InChIKey}')
        if self.PubChem:
            keys.append(f'PubChem={self.PubChem}')
        if self.SMILES:
            keys.append(f'SMILES={self.SMILES}')
        for key in (self.name, self.ID):
            if key and key not in keys:
                keys.append(key)
        return tuple(keys)

    def known_formula(self) -> Optional[str]:
        """Return the formula, read from the InChI if not given."""
        if self.formula:
            return self.formula
        if self.InChI:
            match = _FORMULA_FROM_INCHI.match(self.InChI)
            if match:
                return match.group(1).replace('.', '')
        return None

    def merged(self, **updates) -> ChemicalIdentity:
        """Return a copy where missing fields are filled from `updates`."""
        missing = {
            name: value for name, value in updates.items()
            if value is not None and getattr(self, name) is None
        }
        return replace(self, **missing)

    def to_dict(self) -> dict:
        """Return the fields that are set as a dictionary."""
        return {
            name: value for name, value in asdict(self).items()
            if value is not None and not (name == 'pure' and value)
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ChemicalIdentity:
        """Create an identity from a dictionary, ignoring unknown keys."""
        names = cls.__dataclass_fields__
        return cls(**{i: j for i, j in data.items() if i in names})

    @classmethod
    def from_chemical(cls, chemical: tmo.Chemical) -> ChemicalIdentity:
        """
        Create an identity from a thermosteam chemical.

        Examples
        --------
        >>> import thermosteam as tmo
        >>> from biosteam.interop.chemicals import ChemicalIdentity
        >>> identity = ChemicalIdentity.from_chemical(tmo.Chemical('Water'))
        >>> identity.CAS, identity.InChIKey, identity.pure
        ('7732-18-5', 'XLYOFNOQVPJJNP-UHFFFAOYSA-N', True)

        """
        CAS = chemical.CAS if is_cas_number(chemical.CAS) else None
        iupac_name = chemical.iupac_name
        if isinstance(iupac_name, (tuple, list)):
            iupac_name = iupac_name[0] if iupac_name else None
        pubchem = chemical.pubchemid
        chebi = None
        for alias in sorted(chemical.aliases):
            chebi = (
                normalize_chebi(alias) if 'CHEBI' in alias.upper() else None
            )
            if chebi:
                break
        return cls(
            ID=chemical.ID,
            name=chemical.common_name or None,
            CAS=CAS,
            InChIKey=chemical.InChI_key or None,
            InChI=chemical.InChI or None,
            SMILES=chemical.smiles or None,
            ChEBI=chebi,
            PubChem=str(pubchem) if pubchem else None,
            IUPAC_name=iupac_name or None,
            formula=chemical.formula or None,
            MW=chemical.MW,
            phase=chemical.locked_state,
            pure=CAS is not None,
        )


class ChemicalResolver(Protocol):
    """Protocol of objects that resolve identities to chemicals."""

    def resolve(self, identity: ChemicalIdentity) -> Optional[tmo.Chemical]:
        """Return a chemical for `identity`, or None if unknown."""


def _renamed(chemical: tmo.Chemical, ID: str) -> tmo.Chemical:
    """Return a copy of a chemical with another ID and the same CAS number."""
    return chemical.copy(ID, CAS=chemical.CAS)


def _tag(chemical: tmo.Chemical, identity: ChemicalIdentity) -> tmo.Chemical:
    """Record identifiers thermosteam has no field for as aliases."""
    if identity.ChEBI:
        chemical.aliases.add(identity.ChEBI)
    return chemical


def _complete(chemical: tmo.Chemical, identity: ChemicalIdentity,
              warn: bool = True) -> tmo.Chemical:
    """
    Make a chemical usable in mass and energy balances by fixing its phase
    when it has no vapor pressure model and defaulting missing properties.
    """
    missing = set(chemical.get_missing_properties())
    if not missing & {'V', 'Cn', 'mu', 'kappa', 'Psat', 'Hvap', 'Tb'}:
        return chemical
    if chemical.locked_state is None and missing & {'Psat', 'Tb', 'Hvap'}:
        chemical.at_state(identity.phase or 'l')
        missing = set(chemical.get_missing_properties())
    if 'Cn' in missing and chemical.MW:
        molar_heat_capacity = _WATER_HEAT_CAPACITY * chemical.MW
        if molar_heat_capacity > _MAX_MOLAR_HEAT_CAPACITY:
            # Macromolecules exceed the range of the molar heat capacity
            # models; cap the estimate rather than fail in energy balances.
            Cn = chemical.Cn
            for handle in (getattr(Cn, 'l', None), getattr(Cn, 's', None), Cn):
                if handle is not None and hasattr(handle, 'add_method'):
                    handle.add_method(0.9 * _MAX_MOLAR_HEAT_CAPACITY)
    defaulted = sorted(missing & {'V', 'Cn', 'mu', 'kappa', 'sigma', 'Hf'})
    chemical.default()
    if warn and defaulted:
        warnings.warn(
            f"properties {defaulted} of chemical '{chemical.ID}' are not "
            "available and default to those of water",
            PlaceholderWarning, stacklevel=4,
        )
    return chemical


class UserChemicalsResolver:
    """
    Resolve identities against chemicals the user already defined.

    Chemicals are matched by ID, then by CAS number, then by InChIKey.

    Parameters
    ----------
    chemicals :
        Chemicals to draw from.

    Examples
    --------
    >>> import thermosteam as tmo
    >>> from biosteam.interop.chemicals import (
    ...     ChemicalIdentity, UserChemicalsResolver)
    >>> resolver = UserChemicalsResolver([tmo.Chemical('Water')])
    >>> resolver.resolve(ChemicalIdentity('H2O', CAS='7732-18-5')).ID
    'H2O'
    >>> resolver.resolve(ChemicalIdentity('Ethanol')) is None
    True

    """

    def __init__(self, chemicals: Iterable[tmo.Chemical]):
        self._chemicals = list(chemicals)

    def resolve(self, identity: ChemicalIdentity) -> Optional[tmo.Chemical]:
        """Return the matching user chemical, or None."""
        for chemical in self._chemicals:
            if chemical.ID == identity.ID:
                return chemical
        for chemical in self._chemicals:
            if identity.CAS and chemical.CAS == identity.CAS:
                return _renamed(chemical, identity.ID)
            if identity.InChIKey and chemical.InChI_key == identity.InChIKey:
                return _renamed(chemical, identity.ID)
        return None


class PropertyDatabaseResolver:
    """
    Resolve identities from user-provided property records.

    Parameters
    ----------
    records :
        Keyword arguments for :class:`thermosteam.Chemical` (e.g. `MW`,
        `formula`, `phase`, `Hf`, `rho`, `Cp`) by chemical ID, CAS number,
        InChIKey or ChEBI identifier.

    Notes
    -----
    Records describe chemicals that databases do not know. Properties a
    record does not provide default to those of water.

    Examples
    --------
    >>> from biosteam.interop.chemicals import (
    ...     ChemicalIdentity, PropertyDatabaseResolver)
    >>> resolver = PropertyDatabaseResolver(
    ...     {'Enzyme': {'MW': 1200., 'phase': 'l'}})
    >>> enzyme = resolver.resolve(ChemicalIdentity('Enzyme'))
    >>> enzyme.MW, enzyme.locked_state
    (1200.0, 'l')

    """

    def __init__(self, records: Mapping[str, Mapping[str, Any]]):
        self._records = dict(records)

    def resolve(self, identity: ChemicalIdentity) -> Optional[tmo.Chemical]:
        """Return a chemical built from the matching record, or None."""
        for key in (identity.ID, identity.CAS, identity.InChIKey,
                    identity.ChEBI):
            if key and key in self._records:
                break
        else:
            return None
        kwargs = dict(self._records[key])
        kwargs.setdefault('phase', identity.phase or 'l')
        if 'MW' not in kwargs and 'formula' not in kwargs:
            formula = identity.known_formula()
            if formula:
                kwargs['formula'] = formula
            elif identity.MW:
                kwargs['MW'] = identity.MW
        chemical = tmo.Chemical(identity.ID, search_db=False, **kwargs)
        return _tag(_complete(chemical, identity, warn=False), identity)


class ThermosteamResolver:
    """
    Resolve identities with the databases available to thermosteam.

    The search terms of :meth:`ChemicalIdentity.search_keys` are tried in
    order. Chemicals that are found but lack property models are completed
    with default properties and a :class:`PlaceholderWarning`.

    Examples
    --------
    >>> from biosteam.interop.chemicals import (
    ...     ChemicalIdentity, ThermosteamResolver)
    >>> resolver = ThermosteamResolver()
    >>> identity = ChemicalIdentity(
    ...     'Formate', InChIKey='BDAGIHXWWSANSR-UHFFFAOYSA-N')
    >>> chemical = resolver.resolve(identity)
    >>> chemical.ID, chemical.CAS
    ('Formate', '64-18-6')

    """

    def resolve(self, identity: ChemicalIdentity) -> Optional[tmo.Chemical]:
        """Return the chemical found in the databases, or None."""
        if not identity.pure:
            return None
        for key in identity.search_keys():
            try:
                chemical = tmo.Chemical(identity.ID, search_ID=key)
            except Exception:  # Lookup failures raise a variety of errors.
                continue
            return _tag(_complete(chemical, identity), identity)
        return None


class PubChemResolver:
    """
    Resolve identities online with the PubChem PUG REST service.

    PubChem provides the formula, molecular weight and InChIKey of a
    compound. With these, thermosteam's databases are searched again; if
    the chemical is still unknown, a chemical with default properties is
    created from the formula. ChEBI identifiers are resolved through
    PubChem's synonym index.

    Parameters
    ----------
    fetch :
        Function that returns the body of a URL as text. Defaults to an
        HTTPS request with a timeout; inject a stub for testing or to use a
        cache.
    timeout :
        Timeout of the default request [s].

    Notes
    -----
    This resolver sends chemical identifiers to an external service and is
    therefore never used unless explicitly requested.

    """

    URL = (
        'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{namespace}/'
        '{key}/property/MolecularFormula,MolecularWeight,InChIKey,IUPACName'
        '/JSON'
    )

    def __init__(self, fetch: Optional[Callable[[str], str]] = None,
                 timeout: float = 10.):
        self.timeout = timeout
        self._fetch = fetch or self._request
        self._fallback = ThermosteamResolver()

    def _request(self, url: str) -> str:
        with urllib.request.urlopen(url, timeout=self.timeout) as response:
            return response.read().decode('utf-8')

    def _queries(self, identity: ChemicalIdentity):
        if identity.PubChem:
            yield 'cid', identity.PubChem
        if identity.InChIKey:
            yield 'inchikey', identity.InChIKey
        if identity.ChEBI:
            yield 'name', identity.ChEBI
        if is_cas_number(identity.CAS):
            yield 'name', identity.CAS
        if identity.name:
            yield 'name', identity.name

    def lookup(self, identity: ChemicalIdentity) -> Optional[ChemicalIdentity]:
        """Return `identity` completed with PubChem data, or None."""
        for namespace, key in self._queries(identity):
            url = self.URL.format(
                namespace=namespace, key=urllib.parse.quote(str(key), safe=''),
            )
            try:
                properties = json.loads(self._fetch(url))
                record = properties['PropertyTable']['Properties'][0]
            except (OSError, ValueError, KeyError, IndexError, TypeError):
                continue
            try:
                MW = float(record.get('MolecularWeight'))
            except (TypeError, ValueError):
                MW = None
            CID = record.get('CID')
            return identity.merged(
                PubChem=str(CID) if CID else None,
                InChIKey=record.get('InChIKey'),
                formula=record.get('MolecularFormula'),
                IUPAC_name=record.get('IUPACName'),
                MW=MW,
            )
        return None

    def resolve(self, identity: ChemicalIdentity) -> Optional[tmo.Chemical]:
        """Return a chemical based on PubChem data, or None."""
        if not identity.pure:
            return None
        found = self.lookup(identity)
        if found is None:
            return None
        chemical = self._fallback.resolve(found)
        if chemical is None:
            chemical = PlaceholderResolver().resolve(found)
        return chemical


class PlaceholderResolver:
    """
    Create chemicals with default properties for unknown identities.

    The chemical keeps the formula or molecular weight of the identity and
    takes all other properties from water, so that mass balances are
    meaningful while energy balances are only indicative. A
    :class:`PlaceholderWarning` is issued for each placeholder.

    Examples
    --------
    >>> import warnings
    >>> from biosteam.interop.chemicals import (
    ...     ChemicalIdentity, PlaceholderResolver)
    >>> with warnings.catch_warnings():
    ...     warnings.simplefilter('ignore')
    ...     chemical = PlaceholderResolver().resolve(
    ...         ChemicalIdentity('NADH', formula='C21H29N7O14P2'))
    >>> round(chemical.MW, 2)
    665.44

    """

    def resolve(self, identity: ChemicalIdentity) -> tmo.Chemical:
        """Return a placeholder chemical for `identity`."""
        kwargs = {'phase': identity.phase or 'l'}
        formula = identity.known_formula()
        if formula:
            kwargs['formula'] = formula
        elif identity.MW:
            kwargs['MW'] = identity.MW
        try:
            chemical = tmo.Chemical(identity.ID, search_db=False, **kwargs)
        except Exception:
            # The formula could not be interpreted; keep the weight only.
            kwargs.pop('formula', None)
            if identity.MW:
                kwargs['MW'] = identity.MW
            chemical = tmo.Chemical(identity.ID, search_db=False, **kwargs)
        warnings.warn(
            f"chemical '{identity.ID}' was not found in any database; "
            "using default (water-like) properties",
            PlaceholderWarning, stacklevel=3,
        )
        return _tag(_complete(chemical, identity, warn=False), identity)


@dataclass
class ResolverChain:
    """
    Try resolvers in order until one provides a chemical.

    Parameters
    ----------
    resolvers :
        Resolvers, in order of precedence.

    Examples
    --------
    >>> from biosteam.interop.chemicals import (
    ...     ChemicalIdentity, PropertyDatabaseResolver, ResolverChain,
    ...     ThermosteamResolver)
    >>> chain = ResolverChain([
    ...     PropertyDatabaseResolver({'Enzyme': {'MW': 1200.}}),
    ...     ThermosteamResolver(),
    ... ])
    >>> chemicals = chain.create_chemicals([
    ...     ChemicalIdentity('Water', CAS='7732-18-5'),
    ...     ChemicalIdentity('Enzyme', pure=False),
    ... ])
    >>> chemicals.IDs
    ('Water', 'Enzyme')

    """

    resolvers: list = field(default_factory=list)

    def resolve(self, identity: ChemicalIdentity) -> tmo.Chemical:
        """
        Return a chemical for `identity`.

        Raises
        ------
        ChemicalResolutionError
            If no resolver knows the chemical.

        """
        for resolver in self.resolvers:
            chemical = resolver.resolve(identity)
            if chemical is not None:
                if chemical.ID != identity.ID:
                    chemical = _renamed(chemical, identity.ID)
                return chemical
        raise ChemicalResolutionError(
            f"chemical '{identity.ID}' could not be resolved; provide it "
            "through the 'chemicals' or 'properties' arguments, or allow "
            "placeholders"
        )

    def create_chemicals(
            self, identities: Iterable[ChemicalIdentity]) -> tmo.Chemicals:
        """Return compiled chemicals for all identities, in order."""
        chemicals = tmo.Chemicals([self.resolve(i) for i in identities])
        chemicals.compile()
        return chemicals


def default_resolver(
    chemicals: Optional[Iterable[tmo.Chemical]] = None,
    properties: Optional[Mapping[str, Mapping[str, Any]]] = None,
    online: bool = False,
    placeholders: bool = True,
) -> ResolverChain:
    """
    Return the standard resolver chain.

    Parameters
    ----------
    chemicals :
        Chemicals defined by the user; these take precedence.
    properties :
        Property records for chemicals unknown to databases; see
        :class:`PropertyDatabaseResolver`.
    online :
        Whether to query PubChem for chemicals that thermosteam does not
        know. Identifiers are sent to an external service.
    placeholders :
        Whether to create chemicals with default properties as a last
        resort instead of raising an error.

    """
    resolvers = []
    if chemicals is not None:
        resolvers.append(UserChemicalsResolver(chemicals))
    if properties:
        resolvers.append(PropertyDatabaseResolver(properties))
    resolvers.append(ThermosteamResolver())
    if online:
        resolvers.append(PubChemResolver())
    if placeholders:
        resolvers.append(PlaceholderResolver())
    return ResolverChain(resolvers)
