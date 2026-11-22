"""
Convert the data of Behr et al. (2024) into standard exchange formats.

The script reads the two spreadsheets of the repository
https://github.com/TUDoAD/DWSIM-EnzymeML-KG (CC BY 4.0,
https://doi.org/10.5281/zenodo.10613473) and writes

* ``laccase_abts.enzymeml.json``: the EnzymeML spreadsheet in the
  EnzymeML version 2 JSON serialisation, and
* ``laccase_abts.dexpi.xml``: the process flow diagram of the process
  spreadsheet as a DEXPI 2.0 process model.

Usage::

    git clone https://github.com/TUDoAD/DWSIM-EnzymeML-KG
    python build_inputs.py DWSIM-EnzymeML-KG

Nothing is added to the data. Where the two formats ask for something the
spreadsheets state differently, the choice is listed in README.md.
"""
import json
import os
import sys
import warnings

import openpyxl

from biosteam.interop.dexpi import model as m
from biosteam.interop.dexpi import validation, xmlio

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE = 'https://github.com/TUDoAD/DWSIM-EnzymeML-KG'
SOURCE_DOI = 'https://doi.org/10.5281/zenodo.10613473'
PAPER_DOI = 'https://doi.org/10.3390/pr12030597'


# %% Reading the spreadsheets

def rows(workbook, sheet):
    """Return the non-empty rows of a sheet as lists of cell values."""
    return [
        list(row) for row in workbook[sheet].iter_rows(values_only=True)
        if any(i not in (None, '') for i in row)
    ]


def table(workbook, sheet, key_column=0):
    """Return a sheet with one property per row as a dictionary."""
    result = {}
    for row in rows(workbook, sheet):
        if row[key_column] is not None:
            # Some labels occur twice; the first occurrence is the table.
            result.setdefault(str(row[key_column]).strip(),
                              row[key_column + 1:])
    return result


def text(value):
    return str(value).strip() if value not in (None, '') else None


# %% EnzymeML

def unit(name, *base_units):
    return {'id': name, 'name': name, 'base_units': [
        {'kind': kind, 'exponent': exponent, 'multiplier': 1.0,
         'scale': float(scale)}
        for kind, exponent, scale in base_units
    ]}


UNITS = {
    'ml': unit('ml', ('litre', 1, -3)),
    'mmole / l': unit('mmol / l', ('mole', 1, -3), ('litre', -1, 0)),
    'umole / l': unit('umol / l', ('mole', 1, -6), ('litre', -1, 0)),
    '1 / s': unit('1 / s', ('second', -1, 0)),
    'sec': unit('s', ('second', 1, 0)),
    'K': unit('K', ('kelvin', 1, 0)),
}


def build_enzymeml(workbook, process):
    general = table(workbook, 'General Information')
    vessels = [
        {'id': row[0], 'name': text(row[1]), 'volume': float(row[2]),
         'unit': UNITS[row[3]], 'constant': True}
        for row in rows(workbook, 'Vessels')[2:] if text(row[1])
    ]
    vessel_ids = {i['name']: i['id'] for i in vessels}

    small_molecules = []
    for row in rows(workbook, 'Reactants')[2:]:
        if not text(row[1]):
            continue
        molecule = {
            'id': row[0], 'name': text(row[1]),
            'vessel_id': vessel_ids[text(row[2])],
            'constant': text(row[3]) == 'Constant',
            'synonymous_names': [], 'references': [],
        }
        if text(row[4]):
            molecule['canonical_smiles'] = text(row[4])
        if len(row) > 5 and text(row[5]):
            molecule['inchi'] = text(row[5])
        small_molecules.append(molecule)
    # Water takes part in the reaction according to the process
    # spreadsheet but is not a species of the EnzymeML spreadsheet.
    small_molecules.append({
        'id': 's6', 'name': 'water', 'vessel_id': vessels[0]['id'],
        'constant': True, 'inchi': 'InChI=1S/H2O/h1H2',
        'canonical_smiles': 'O', 'synonymous_names': [], 'references': [],
    })

    proteins = []
    for row in rows(workbook, 'Proteins')[2:]:
        if not text(row[1]):
            continue
        proteins.append({
            'id': row[0], 'name': text(row[1]), 'sequence': text(row[2]),
            'vessel_id': vessel_ids[text(row[3])],
            'constant': text(row[4]) == 'Constant',
            'organism': text(row[5]), 'ecnumber': text(row[6]),
            'references': [
                f'https://www.uniprot.org/uniprotkb/{text(row[7])}'],
        })

    def in_vessel(items, name, vessel_id):
        for item in items:
            if item['name'] == name and item['vessel_id'] == vessel_id:
                return item['id']
        raise KeyError((name, vessel_id))

    # Stoichiometric numbers come from the process spreadsheet
    # (4 ABTS_red + O2 -> 4 ABTS_ox + 2 H2O) and are written per mole of
    # ABTS, because the rate law gives the rate of ABTS oxidation.
    coefficients = dict(zip(
        process['substances']['Name'],
        process['substances']['hasStoichiometricCoefficient'],
    ))
    per_abts = abs(coefficients['ABTS_red'])
    numbers = {
        'ABTS reduced': abs(coefficients['ABTS_red']) / per_abts,
        'oxygen': abs(coefficients['Oxygen']) / per_abts,
        'ABTS oxidized': abs(coefficients['ABTS_ox']) / per_abts,
        'water': abs(coefficients['Water']) / per_abts,
    }

    equations = {
        text(row[0]): text(row[1])
        for row in rows(workbook, 'Kinetic Models')[2:]
        if text(row[0]) and text(row[1]) and '/' in str(row[1])
    }
    equations = {i.lower(): j for i, j in equations.items()}
    reactions = []
    conditions = {}
    for row in rows(workbook, 'Reactions')[2:]:
        if not text(row[1]):
            continue
        name = text(row[1])
        protein = next(i for i in proteins if i['name'] == text(row[8])
                       or name.endswith(i['name'].split(', ')[-1]))
        vessel_id = protein['vessel_id']
        temperature = float(row[2]) + 273.15  # degC in the spreadsheet
        conditions[name] = (temperature, float(row[4]))
        reactions.append({
            'id': row[0], 'name': name,
            'reversible': text(row[5]) != 'irreversible',
            'kinetic_law': {
                'species_id': f'v_{row[0]}',
                'equation': equations[name.lower().replace('-', ' ')],
                'equation_type': 'rateLaw', 'variables': [],
            },
            'reactants': [
                {'species_id': in_vessel(small_molecules, i, vessel_id),
                 'stoichiometry': numbers[i]}
                for i in ('ABTS reduced', 'oxygen')
            ],
            'products': [
                {'species_id': in_vessel(small_molecules, 'ABTS oxidized',
                                         vessel_id),
                 'stoichiometry': numbers['ABTS oxidized']},
                {'species_id': 's6', 'stoichiometry': numbers['water']},
            ],
            'modifiers': [
                {'species_id': protein['id'], 'role': 'biocatalyst'},
            ],
        })

    parameters = []
    for row in rows(workbook, 'Kinetic Models')[2:]:
        if len(row) > 6 and text(row[4]) in ('Km', 'kcat'):
            parameters.append({
                'id': text(row[4]), 'name': text(row[4]),
                'symbol': text(row[4]), 'value': float(row[5]),
                'unit': UNITS[text(row[6])], 'fit': False, 'constant': True,
            })

    measurements = build_measurements(
        workbook, reactions, proteins, small_molecules, conditions)
    return {
        'name': text(general['Title'][0]),
        'version': '2',
        'description': (
            'Laccase-catalysed oxidation of ABTS in a helically coiled and '
            'a straight capillary flow reactor. Rewritten in the EnzymeML '
            f'version 2 JSON serialisation from the EnzymeML spreadsheet of '
            f'Behr et al., Processes 12 (2024) 597, {PAPER_DOI}; data from '
            f'{SOURCE} ({SOURCE_DOI}), licensed under CC BY 4.0.'
        ),
        'created': str(general['Date of creation'][0])[:10],
        'creators': [],
        'vessels': vessels,
        'proteins': proteins,
        'complexes': [],
        'small_molecules': small_molecules,
        'reactions': reactions,
        'measurements': measurements,
        'equations': [],
        'parameters': parameters,
        'references': [PAPER_DOI, SOURCE_DOI],
    }


def build_measurements(workbook, reactions, proteins, small_molecules,
                       conditions):
    """Return one measurement per experiment of the 'Data' sheet."""
    reaction_by_name = {i['name']: i for i in reactions}
    protein_by_name = {i['name']: i for i in proteins}
    experiments = {}
    for row in rows(workbook, 'Data'):
        if not isinstance(row[0], (int, float)):
            continue
        experiment = experiments.setdefault(int(row[0]), {
            'reaction': text(row[2]), 'protein': text(row[3]), 'species': [],
        })
        series = [float(i) for i in row[9:] if isinstance(i, (int, float))]
        if str(row[1]).startswith('Time'):
            experiment['time'] = series
            experiment['time_unit'] = UNITS[text(row[8])]
        else:
            experiment['enzyme'] = (float(row[4]), UNITS[text(row[5])])
            experiment['species'].append(
                (text(row[6]), float(row[7]), UNITS[text(row[8])], series))
    measurements = []
    for number, experiment in sorted(experiments.items()):
        protein = protein_by_name[experiment['protein']]
        vessel_id = protein['vessel_id']
        temperature, pH = conditions[experiment['reaction']]
        value, enzyme_unit = experiment['enzyme']
        species_data = [{
            'species_id': protein['id'], 'initial': value,
            'data_unit': enzyme_unit, 'data': [], 'time': [],
            'data_type': 'concentration', 'is_simulated': False,
        }]
        for name, initial, data_unit, series in experiment['species']:
            species_id = next(
                i['id'] for i in small_molecules
                if i['name'] == name and i['vessel_id'] == vessel_id)
            entry = {
                'species_id': species_id, 'initial': initial,
                'data_unit': data_unit, 'data': series,
                'time': experiment['time'] if series else [],
                'data_type': 'concentration', 'is_simulated': False,
            }
            if series:
                entry['time_unit'] = experiment['time_unit']
            species_data.append(entry)
        measurements.append({
            'id': f'm{number}',
            'name': f"Experiment {number}: {experiment['reaction']}",
            'species_data': species_data,
            'group_id': reaction_by_name[experiment['reaction']]['id'],
            'ph': pH, 'temperature': temperature,
            'temperature_unit': UNITS['K'],
        })
    return measurements


# %% DEXPI

def quantity(value, enumeration, literal, identifier):
    obj = m.qualified_value(m.physical_quantity(value, enumeration, literal))
    obj.id = identifier
    return obj


def build_dexpi(process):
    substances = process['substances']
    streams = process['streams']
    reactor = process['reactor']
    names = [text(i) for i in streams['Name'] if text(i)]
    properties = process['properties']
    cas = dict(zip(properties['Name'], properties['hasCAS_Number']))
    in_database = dict(zip(substances['Name'], substances['inDWSIMdatabase']))

    root = m.DexpiObject(m.ENGINEERING_MODEL, 'model')
    model = root.add('ConceptualModel',
                     m.DexpiObject(m.PROCESS_MODEL, 'process'))
    components = []
    for name in names:
        role = dict(zip(substances['Name'], substances['has_role']))[name]
        kind = ('CustomMaterialComponent' if role == 'catalyst'
                else 'PureMaterialComponent')
        component = m.DexpiObject(
            m.process_type(kind), f'component-{name}',
            {'Identifier': name, 'Label': name})
        if name == 'Water':
            cas[name] = '7732-18-5'
        elif name == 'Oxygen':
            cas[name] = '7782-44-7'
        if cas.get(name) and (in_database[name] or kind.startswith('Pure')):
            component.add('PersistentIdentifiers', m.DexpiObject(
                m.PERSISTENT_IDENTIFIER, f'cas-{name}',
                {'Context': 'CAS', 'Value': str(cas[name])}))
        components.append(model.add('MaterialComponents', component))
    component_list = model.add('ListsOfMaterialComponents', m.DexpiObject(
        m.process_type('ListOfMaterialComponents'), 'components'))
    component_list.refer('Component', *components)
    template = model.add('MaterialTemplates', m.DexpiObject(
        m.process_type('MaterialTemplate'), 'material', {
            'Description': m.multi_language_string(
                'Gas-liquid reaction mixture'),
            'Identifier': 'mixture', 'Label': 'mixture',
            'NumberOfMaterialComponents': len(names),
            'NumberOfPhases': 2, 'PhaseLabel': ['l', 'g'],
        }))
    template.refer('ListOfComponents', component_list)

    flows = dict(zip(names, streams['hasCompoundMolarFlow']))  # mol/s
    enters = dict(zip(names, streams['EntersAtObject']))
    temperature = float(streams['hasTemperature'][0])  # degC

    objects, connections = {}, {}
    for row in process['pfd'][1:]:
        name = text(row[0])
        objects[name] = {
            'type': text(row[1]), 'to': text(row[3]),
            'volume_flow': row[4] if len(row) > 4 else None,
        }

    def add_stream(name):
        info = objects[name]
        stream = model.add('ProcessConnections', m.DexpiObject(
            m.process_type('Stream'), f'stream-{name}',
            {'Identifier': name, 'Label': name}))
        stream.refer('MaterialTemplateReference', template)
        if info['volume_flow']:
            # Feed stream: conditions and composition are specified.
            phase = 'g' if enters.get('Oxygen') == name else 'l'
            stream.add('Temperature', quantity(
                temperature, 'TemperatureUnit', 'DegreeCelsius',
                f'{name}-temperature'))
            stream.add('Pressure', quantity(
                float(reactor['hasPressureValue'][0]), 'PressureAbsoluteUnit',
                'Pascal', f'{name}-pressure'))
            stream.add('VolumeFlow', quantity(
                float(info['volume_flow']), 'VolumeFlowRateUnit',
                'MetreCubedPerSecond', f'{name}-volumeflow'))
            state = model.add('MaterialStates', m.DexpiObject(
                m.process_type('MaterialState'), f'{name}-state', {
                    'Description': f'Composition of {name}',
                    'Identifier': f'{name}.state', 'Label': name}))
            total = state.add('State', m.DexpiObject(
                m.process_type('MaterialStateType'), f'{name}-total', {
                    'Description': m.multi_language_string('Total'),
                    'Identifier': f'{name}.total', 'Label': 'total'}))
            composition = total.add('Composition', m.DexpiObject(
                m.process_type('Composition'), f'{name}-composition', {
                    'Display': m.process_enumeration(
                        'CompositionDisplay', 'AbsoluteValue')}))
            for i, component in enumerate(names):
                flow = (float(flows[component])
                        if enters[component] == name else 0.)
                composition.add('MoleFlow', quantity(
                    flow / 1000., 'MoleFlowRateUnit', 'KilomolePerSecond',
                    f'{name}-moleflow-{i}'))
            state.add('Phases', m.DexpiObject(
                m.process_type('MaterialStateType'), f'{name}-phase', {
                    'Description': m.multi_language_string('Phase'),
                    'Identifier': f'{name}.{phase}', 'Label': phase}))
            stream.refer('MaterialStateReference', state)
        connections[name] = stream
        return stream

    def add_step(cls, identifier, label, ports, **data):
        step = model.add('ProcessSteps', m.DexpiObject(
            m.process_type(cls), f'step-{identifier}',
            {'Identifier': identifier, 'Label': label, **data}))
        for port_name, direction, connection in ports:
            port = step.add('Ports', m.DexpiObject(
                m.process_type('MaterialPort'), f'{identifier}-{port_name}', {
                    'Identifier': port_name,
                    'NominalDirection': m.process_enumeration(
                        'PortDirection', direction)}))
            port.refer('ConnectorReference', connection)
            end = 'Source' if direction == 'Outlet' else 'Target'
            connection.refer(end, port)
        return step

    material = [i for i, j in objects.items() if j['type'] == 'MaterialStream']
    for name in material:
        add_stream(name)
    for name in material:  # Feeds and products end at the battery limit
        if not any(j['to'] == name for j in objects.values()):
            add_step('Source', f'source_{name}', name,
                     [('out0', 'Outlet', connections[name])])
        if not objects[name]['to']:
            add_step('Sink', f'sink_{name}', name,
                     [('in0', 'Inlet', connections[name])])
    for name, info in objects.items():
        if info['type'] not in ('Mixer', 'RCT_PFR'):
            continue  # Streams; the energy stream is not represented
        inlets = [i for i in material if objects[i]['to'] == name]
        ports = [(f'in{i}', 'Inlet', connections[j])
                 for i, j in enumerate(inlets)]
        ports.append(('out0', 'Outlet', connections[info['to']]))
        if info['type'] == 'Mixer':
            add_step('MixingSimple', name, 'Gas-liquid mixing', ports)
            continue
        step = add_step(
            'ReactingChemicals', name, 'Capillary flow reactor', ports,
            Method=m.process_enumeration('ReactionProcessType', 'Tubular'))
        step.add('Temperature', quantity(
            temperature, 'TemperatureUnit', 'DegreeCelsius',
            f'{name}-temperature'))
        step.add('Pressure', quantity(
            float(reactor['hasPressureValue'][0]), 'PressureAbsoluteUnit',
            'Pascal', f'{name}-pressure'))
        step.add('Volume', quantity(
            float(reactor['hasVolumeValue'][0]), 'VolumeUnit', 'MetreCubed',
            f'{name}-volume'))
    root.data['OriginatingSystemName'] = (
        'Converted from the process spreadsheet of Behr et al. (2024)')
    document = m.DexpiDocument(
        root, name='laccase_abts_flow_reactor', uri=SOURCE_DOI)
    validation.validate(document)
    return document


# %% Main

def main(repository):
    warnings.simplefilter('ignore')
    folder = os.path.join(repository, 'ELNs')
    enzymeml_workbook = openpyxl.load_workbook(
        os.path.join(folder, 'EnzymeML_Template_18-8-2021_KR.xlsm'),
        data_only=True, read_only=True, keep_vba=False)
    process_workbook = openpyxl.load_workbook(
        os.path.join(folder, 'New-ELN_Kinetik_1.xlsx'),
        data_only=True, read_only=True)
    process = {
        'substances': table(process_workbook, 'Substances and Reactions'),
        'properties': table(process_workbook, 'Properties for JSON-file'),
        'streams': table(process_workbook, 'Material Streams'),
        'reactor': table(process_workbook, 'Reactor Specification'),
        'pfd': rows(process_workbook, 'PFD'),
    }
    enzymeml = build_enzymeml(enzymeml_workbook, process)
    path = os.path.join(HERE, 'laccase_abts.enzymeml.json')
    with open(path, 'w', encoding='utf-8') as file:
        json.dump(enzymeml, file, indent=2)
        file.write('\n')
    xmlio.write(build_dexpi(process),
                os.path.join(HERE, 'laccase_abts.dexpi.xml'))
    print(f"{len(enzymeml['measurements'])} measurements, "
          f"{len(enzymeml['reactions'])} reactions written")


if __name__ == '__main__':
    main(sys.argv[1])
