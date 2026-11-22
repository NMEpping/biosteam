"""
First use case: NADH regeneration by formate dehydrogenase in a stirred
tank, simulated in BioSTEAM from an EnzymeML document and a DEXPI
flowsheet.

Run ``python run.py`` in this folder. The script prints the results
reported in the manuscript and writes the simulated system as
``fdh_simulated.dexpi.xml``.
"""
import os
import warnings

import biosteam as bst
import numpy as np
from biosteam.interop import chemicals, dexpi, enzymeml

HERE = os.path.dirname(os.path.abspath(__file__))
ENZYMEML = os.path.join(HERE, 'fdh_nadh_regeneration.enzymeml.json')
DEXPI = os.path.join(HERE, 'fdh_flowsheet.dexpi.xml')

#: What chemical databases do not provide: the formula of NADH and the
#: approximate molecular weight [g/mol] of the enzyme.
PROPERTIES = {
    'nadh': {'formula': 'C21H29N7O14P2'},
    'fdh': {'MW': 40400.},
}


def main():
    warnings.simplefilter('ignore', chemicals.PlaceholderWarning)
    warnings.simplefilter('ignore', dexpi.DexpiWarning)
    warnings.simplefilter('ignore', bst.exceptions.CostWarning)
    document = enzymeml.read(ENZYMEML)
    bst.main_flowsheet.set_flowsheet('fdh')
    bst.settings.set_thermo(document.create_chemicals(properties=PROPERTIES))
    system = bst.System.from_dexpi(DEXPI)
    reactor = system.flowsheet.unit.R1
    document.parameterize(reactor)
    reactor.tau = 2.  # Residence time [hr]
    system.simulate()
    feed = system.flowsheet.stream.feed
    product = system.flowsheet.stream.product
    conversion = 1 - product.imol['nad'] / feed.imol['nad']
    print(f"units: {[f'{i.ID} ({type(i).__name__})' for i in system.units]}")
    print(f"flow pattern from DEXPI: {reactor.flow_pattern}")
    print(f"NAD+ conversion: {conversion:.1%}")
    print(f"NADH produced: {product.imol['nadh'] * 1e3:.2f} mol/hr")
    print(f"reactor volume: "
          f"{reactor.design_results['Reactor volume']:.2f} m3")

    path = os.path.join(HERE, 'fdh_simulated.dexpi.xml')
    system.to_dexpi(path)
    copy = bst.System.from_dexpi(
        path, flowsheet='fdh_copy',
        resolver=chemicals.default_resolver(properties=PROPERTIES))
    copy.simulate()
    same = np.allclose(copy.flowsheet.stream.product.mol, product.mol,
                       rtol=1e-9, atol=1e-18)
    print(f"round trip reproduces the product stream: {bool(same)}")


if __name__ == '__main__':
    main()
