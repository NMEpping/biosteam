"""
Second use case: the laccase/ABTS flow reactor of Behr et al. (2024),
simulated in BioSTEAM from an EnzymeML document and a DEXPI flowsheet.

Run ``python run.py`` in this folder after ``build_inputs.py``. The script
writes ``results.csv`` and ``fig_laccase.pdf`` and prints a summary.
"""
import csv
import json
import os
import warnings

import biosteam as bst
import numpy as np
from biosteam.interop import chemicals, dexpi, enzymeml
from biosteam.interop.dexpi import quantities, xmlio

HERE = os.path.dirname(os.path.abspath(__file__))
ENZYMEML = os.path.join(HERE, 'laccase_abts.enzymeml.json')
DEXPI = os.path.join(HERE, 'laccase_abts.dexpi.xml')

#: Chemical of the flowsheet for each species of the EnzymeML document.
#: The document lists every compound once per vessel; the flowsheet names
#: each compound once.
CHEMICAL_IDS = {
    's0': 'ABTS_red', 's3': 'ABTS_red',
    's1': 'Oxygen', 's4': 'Oxygen',
    's2': 'ABTS_ox', 's5': 'ABTS_ox',
    's6': 'Water',
    'p0': 'Laccase', 'p1': 'Laccase', 'p2': 'Laccase',
}

#: Molecular weights [g/mol] as stated in the process spreadsheet of
#: Behr et al.; ABTS and the enzyme are not in the property databases.
PROPERTIES = {
    'ABTS_red': {'MW': 514.619},
    'ABTS_ox': {'MW': 513.619},
    'Laccase': {'MW': 56000.},
}

#: Experiments at the conditions of the process spreadsheet (11.6 uM
#: laccase, 0.97 mM ABTS, 1.34 mM oxygen), by reactor.
EXPERIMENTS = {'m3': 'helically coiled capillary', 'm9': 'straight capillary'}

#: Concentration of oxidised ABTS in the liquid product [mol/m3] that Behr
#: et al. obtained with DWSIM, as stored in their knowledge graph
#: (ontologies/KG-DWSIM_EnzML_ELN_output.owl, indv_Product_1_ABTS_ox).
BEHR_DWSIM = 0.44847


def load_system():
    """Return the system, its reactor and the EnzymeML document."""
    document = enzymeml.read(ENZYMEML).with_chemical_IDs(CHEMICAL_IDS)
    bst.main_flowsheet.set_flowsheet('laccase')
    bst.settings.set_thermo(document.create_chemicals(properties=PROPERTIES))
    system = bst.System.from_dexpi(DEXPI, chemicals=bst.settings.get_thermo())
    reactor = system.flowsheet.unit.Reactor
    # The reactor is tubular in the DEXPI file, so plug flow is the default.
    document.parameterize(reactor, 'r2')
    return system, reactor, document


def stated_residence_time():
    """Return the residence time [s] from the reactor volume and the
    volumetric feed flows stated in the DEXPI file."""
    process = xmlio.read(DEXPI).process_model
    flow = sum(
        quantities.first(i, 'VolumeFlow', 'volume flow') or 0.
        for i in process.children('ProcessConnections')
    )  # m3/hr
    reactor = next(i for i in process.children('ProcessSteps')
                   if i.get('Identifier') == 'Reactor')
    return quantities.first(reactor, 'Volume', 'volume') / flow * 3600.


def product_concentration(system, reactor, seconds):
    """Return the ABTS_ox concentration [mol/m3] at a residence time."""
    reactor.tau = seconds / 3600.
    system.run()
    product = reactor.outs[0]
    return product.imol['ABTS_ox'] / product.F_vol * 1e3


def measured():
    """Return the measured time courses {id: (time, concentration)}."""
    with open(ENZYMEML, encoding='utf-8') as file:
        data = json.load(file)
    courses = {}
    for measurement in data['measurements']:
        if measurement['id'] not in EXPERIMENTS:
            continue
        series = next(i for i in measurement['species_data'] if i['data'])
        courses[measurement['id']] = (
            np.array(series['time']), np.array(series['data']))
    return courses


def main():
    warnings.simplefilter('ignore', chemicals.PlaceholderWarning)
    warnings.simplefilter('ignore', dexpi.DexpiWarning)
    warnings.simplefilter('ignore', bst.exceptions.CostWarning)
    system, reactor, document = load_system()
    tau = stated_residence_time()
    feed = reactor.ins[0]
    # Mass and energy balances only: BioSTEAM's vessel design and costing
    # for stirred tanks do not apply to a capillary of a few millilitres.
    system.run()
    S0 = feed.imol['ABTS_red'] / feed.F_vol * 1e3
    E0 = feed.imol['Laccase'] / feed.F_vol * 1e3
    print(f"units: {[f'{i.ID} ({type(i).__name__})' for i in system.units]}")
    print(f"flow pattern from DEXPI: {reactor.flow_pattern}")
    print(f"residence time from DEXPI volume and flows: {tau:.1f} s")
    print(f"reactor inlet: {S0:.3f} mM ABTS, {E0 * 1e3:.2f} uM laccase")

    times = np.linspace(0., 62., 63)
    times[0] = 1e-6
    plug = np.array([product_concentration(system, reactor, t) for t in times])

    # Independent check: integrated Michaelis-Menten equation,
    # Km ln(S0/S) + S0 - S = kcat E0 t  (concentrations in mol/m3).
    Km = document.parameters['Km'].value
    kcat = document.parameters['kcat'].value
    S = S0 - plug
    residual = Km * np.log(S0 / S) + S0 - S - kcat * E0 * times
    times[0] = 0.
    print(f"largest deviation from the integrated rate equation: "
          f"{abs(residual).max():.1e} mol/m3")

    document.parameterize(reactor, 'r2', flow_pattern='mixed')
    mixed = np.array(
        [product_concentration(system, reactor, t) for t in times])
    document.parameterize(reactor, 'r2', flow_pattern='plug')
    at_tau = product_concentration(system, reactor, tau)
    print(f"ABTS_ox at {tau:.0f} s: {at_tau:.3f} mM (plug flow)")
    print(f"DWSIM result of Behr et al.: {BEHR_DWSIM:.3f} mM "
          f"({at_tau / BEHR_DWSIM - 1:+.0%} here)")

    courses = measured()
    for identifier, (time, data) in courses.items():
        model = np.interp(time, times, plug)
        early = time < 40.
        print(f"{EXPERIMENTS[identifier]} ({identifier}): measured "
              f"{data[-1]:.3f} mM at {time[-1]:.0f} s; mean absolute "
              f"deviation {abs(model - data)[early].mean():.3f} mM below "
              f"40 s, {abs(model - data)[~early].mean():.3f} mM above")

    # Round trip through DEXPI
    xml = system.to_dexpi()
    copy = bst.System.from_dexpi(xml, flowsheet='laccase_copy')
    copy.run()
    same = np.allclose(copy.flowsheet.stream.Product_1.mol,
                       reactor.outs[0].mol, rtol=1e-9, atol=1e-18)
    print(f"round trip reproduces the product stream: {bool(same)}")

    with open(os.path.join(HERE, 'results.csv'), 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['residence_time_s', 'ABTS_ox_plug_flow_mM',
                         'ABTS_ox_stirred_tank_mM'])
        writer.writerows(zip(times, plug.round(6), mixed.round(6)))
    plot(times, plug, mixed, courses, tau)


def plot(times, plug, mixed, courses, tau):
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    folder = '/usr/share/texmf/fonts/opentype/public/lm'
    if os.path.isdir(folder):  # Latin Modern, as in the manuscript
        for name in ('lmroman10-regular.otf', 'lmroman10-bold.otf'):
            font_manager.fontManager.addfont(os.path.join(folder, name))
        plt.rcParams['font.family'] = 'Latin Modern Roman'
    plt.rcParams.update({
        'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 8,
        'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'mathtext.fontset': 'cm', 'axes.linewidth': 0.6,
        'axes.spines.top': False, 'axes.spines.right': False,
    })
    turquoise, purple, red, gray = '#00C1D4', '#7200FE', '#FF4F4F', '#A8968C'
    fig, ax = plt.subplots(figsize=(390 / 72.27, 2.6))
    ax.plot(times, plug, color=purple, lw=1.6,
            label='BioSTEAM, plug flow (as in DEXPI file)')
    ax.plot(times, mixed, color=gray, lw=1.2, ls='--',
            label='BioSTEAM, stirred tank (for contrast)')
    markers = {'m3': 'o', 'm9': 's'}
    colors = {'m3': turquoise, 'm9': red}
    for identifier, (time, data) in courses.items():
        ax.plot(time, data, markers[identifier], ms=4, mfc='white',
                mec=colors[identifier], mew=1.2,
                label=f'measured, {EXPERIMENTS[identifier]}')
    ax.plot([tau], [BEHR_DWSIM], '*', ms=9, color='#FFDE36', mec='black',
            mew=0.6, label='DWSIM (Behr et al.)')
    ax.axvline(tau, color='black', lw=0.6, ls=':')
    ax.text(tau - 1, 0.02, 'residence time\nfrom DEXPI', ha='right',
            va='bottom', fontsize=8)
    ax.set_xlabel('Residence time [s]')
    ax.set_ylabel('Oxidised ABTS [mM]')
    ax.set_xlim(0, 63)
    ax.set_ylim(0, None)
    ax.legend(frameon=False, loc='upper left', handlelength=2.2)
    fig.tight_layout(pad=0.3)
    fig.savefig(os.path.join(HERE, 'fig_laccase.pdf'))


if __name__ == '__main__':
    main()
