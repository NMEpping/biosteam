# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Tests for the resolution of chemical identities to thermosteam chemicals.
"""
import json
import warnings

import pytest
import thermosteam as tmo
from biosteam.interop.chemicals import (
    ChemicalIdentity,
    ChemicalResolutionError,
    PlaceholderResolver,
    PlaceholderWarning,
    PropertyDatabaseResolver,
    PubChemResolver,
    ResolverChain,
    ThermosteamResolver,
    UserChemicalsResolver,
    default_resolver,
)


def test_identity_search_keys_and_dictionary_round_trip():
    identity = ChemicalIdentity(
        'EtOH', name='ethanol', CAS='64-17-5',
        InChIKey='LFQSCWFLJHTTHZ-UHFFFAOYSA-N', PubChem='702',
    )
    assert identity.search_keys() == (
        '64-17-5', 'InChIKey=LFQSCWFLJHTTHZ-UHFFFAOYSA-N', 'PubChem=702',
        'ethanol', 'EtOH',
    )
    assert ChemicalIdentity.from_dict(identity.to_dict()) == identity
    assert ChemicalIdentity.from_dict({'ID': 'X', 'unknown': 1}) == ChemicalIdentity('X')
    # An identifier that is not a CAS number is not searched as one
    assert ChemicalIdentity('X', CAS='not-a-cas').search_keys() == ('X',)
    assert ChemicalIdentity('X', InChI='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3').known_formula() == 'C2H6O'
    assert identity.merged(name='other', formula='C2H6O') == ChemicalIdentity(
        'EtOH', name='ethanol', CAS='64-17-5', formula='C2H6O',
        InChIKey='LFQSCWFLJHTTHZ-UHFFFAOYSA-N', PubChem='702',
    )


def test_identity_from_chemical():
    identity = ChemicalIdentity.from_chemical(tmo.Chemical('Ethanol'))
    assert identity.CAS == '64-17-5' and identity.pure
    assert identity.InChIKey == 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
    custom = tmo.Chemical('Cells', search_db=False, phase='s', default=True,
                          formula='CH1.8O0.5N0.2')
    identity = ChemicalIdentity.from_chemical(custom)
    assert identity.CAS is None and not identity.pure
    assert identity.phase == 's' and identity.formula


@pytest.mark.parametrize('identity, CAS', [
    (ChemicalIdentity('A', CAS='64-17-5'), '64-17-5'),
    (ChemicalIdentity('A', InChIKey='XLYOFNOQVPJJNP-UHFFFAOYSA-N'), '7732-18-5'),
    (ChemicalIdentity('A', name='acetic acid'), '64-19-7'),
    (ChemicalIdentity('Glucose'), '50-99-7'),
    # The most specific identifier wins over a misleading name
    (ChemicalIdentity('Water', name='ethanol', CAS='7732-18-5'), '7732-18-5'),
])
def test_thermosteam_resolver_fallback_order(identity, CAS):
    chemical = ThermosteamResolver().resolve(identity)
    assert chemical.ID == identity.ID and chemical.CAS == CAS


def test_thermosteam_resolver_skips_unknown_and_custom_components():
    resolver = ThermosteamResolver()
    assert resolver.resolve(ChemicalIdentity('NotAChemicalName123')) is None
    assert resolver.resolve(ChemicalIdentity('Water', pure=False)) is None


def test_incomplete_database_chemicals_are_completed():
    # NAD+ is in the databases, but without models for most properties
    with pytest.warns(PlaceholderWarning, match="'nad'"):
        chemical = ThermosteamResolver().resolve(ChemicalIdentity('nad', CAS='53-84-9'))
    assert chemical.locked_state == 'l'
    assert chemical.formula == 'C21H27N7O14P2'
    stream = tmo.Stream(None, thermo=tmo.Thermo([tmo.Chemical('Water'), chemical]),
                        Water=1., nad=1e-3)
    assert stream.F_vol > 0. and stream.H == pytest.approx(0.)
    stream.T = 310.
    assert stream.H > 0.


def test_user_chemicals_take_precedence():
    mine = tmo.Chemical('Water')
    mine_by_cas = tmo.Chemical('Ethanol')
    resolver = default_resolver(chemicals=[mine, mine_by_cas])
    assert resolver.resolve(ChemicalIdentity('Water')) is mine
    copy = resolver.resolve(ChemicalIdentity('EtOH', CAS='64-17-5'))
    assert copy.ID == 'EtOH' and copy.CAS == '64-17-5'
    assert UserChemicalsResolver([mine]).resolve(ChemicalIdentity('Other')) is None


def test_property_database_and_macromolecules():
    resolver = PropertyDatabaseResolver({
        'CHEBI:16908': {'formula': 'C21H29N7O14P2'},
        'enzyme': {'MW': 40400., 'phase': 'l'},
    })
    nadh = resolver.resolve(ChemicalIdentity('NADH', ChEBI='CHEBI:16908'))
    assert nadh.MW == pytest.approx(665.44, abs=0.01)
    assert 'CHEBI:16908' in nadh.aliases
    assert resolver.resolve(ChemicalIdentity('other')) is None
    # Macromolecules must remain usable in energy balances
    enzyme = resolver.resolve(ChemicalIdentity('enzyme', pure=False))
    stream = tmo.Stream(None, thermo=tmo.Thermo([tmo.Chemical('Water'), enzyme]),
                        Water=50., enzyme=1e-6, T=320.)
    assert stream.H > 0. and stream.F_mass == pytest.approx(50 * 18.01528 + 0.0404)


def test_placeholders_and_strict_resolution():
    identity = ChemicalIdentity('Mystery', MW=250., phase='s')
    with pytest.warns(PlaceholderWarning, match="'Mystery' was not found"):
        chemical = default_resolver().resolve(identity)
    assert chemical.MW == 250. and chemical.locked_state == 's'
    with pytest.raises(ChemicalResolutionError, match="'Mystery' could not be resolved"):
        default_resolver(placeholders=False).resolve(identity)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # A formula that cannot be parsed does not prevent a placeholder
        chemical = PlaceholderResolver().resolve(
            ChemicalIdentity('Odd', formula='not a formula', MW=10.))
    assert chemical.MW == 10.


def test_resolver_chain_creates_compiled_chemicals():
    chain = ResolverChain([
        PropertyDatabaseResolver({'Enzyme': {'MW': 1000.}}),
        ThermosteamResolver(),
    ])
    chemicals = chain.create_chemicals([
        ChemicalIdentity('H2O', CAS='7732-18-5', ChEBI='CHEBI:15377'),
        ChemicalIdentity('Enzyme', pure=False),
    ])
    assert chemicals.IDs == ('H2O', 'Enzyme')
    assert chemicals['CHEBI:15377'] is chemicals.H2O  # ChEBI id as alias


def test_pubchem_resolver_with_stubbed_service():
    requests = []
    record = {'PropertyTable': {'Properties': [{
        'CID': 439153, 'MolecularFormula': 'C21H29N7O14P2',
        'MolecularWeight': '665.4', 'InChIKey': 'TESTKEY-FOR-STUB-N',
        'IUPACName': 'a long name',
    }]}}

    def fetch(url):
        requests.append(url)
        if 'CHEBI%3A16908' not in url:
            raise OSError('not found')
        return json.dumps(record)

    resolver = PubChemResolver(fetch=fetch)
    identity = ChemicalIdentity('NADH', InChIKey='UNKNOWNKEY', ChEBI='CHEBI:16908')
    found = resolver.lookup(identity)
    # Identifiers are tried in order of specificity; the ChEBI id matches
    assert [i.split('/compound/')[1].split('/')[0] for i in requests] == ['inchikey', 'name']
    assert found.formula == 'C21H29N7O14P2' and found.PubChem == '439153'
    assert found.InChIKey == 'UNKNOWNKEY'  # Given identifiers are kept
    with pytest.warns(PlaceholderWarning):
        chemical = resolver.resolve(identity)
    assert chemical.ID == 'NADH' and chemical.MW == pytest.approx(665.44, abs=0.01)
    # Failures of the service fall through to the next resolver
    assert PubChemResolver(fetch=lambda url: 'not json').resolve(identity) is None
    assert resolver.resolve(ChemicalIdentity('enzyme', pure=False)) is None


def test_online_lookup_is_opt_in():
    kinds = [type(i) for i in default_resolver().resolvers]
    assert PubChemResolver not in kinds
    kinds = [type(i) for i in default_resolver(online=True).resolvers]
    assert kinds.index(PubChemResolver) > kinds.index(ThermosteamResolver)
