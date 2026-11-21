# -*- coding: utf-8 -*-
"""
NADH regeneration by formate dehydrogenase in a CSTR: from EnzymeML
kinetics to a DEXPI process model and back.

Run this script with ``python dexpi_fdh_cstr.py``. It needs BioSTEAM only.

"""
import os
import tempfile
import warnings

import biosteam as bst
from biosteam.interop import chemicals, enzymeml, examples


def main(folder=None):
    """Run the workflow and return the original and the re-imported system."""
    warnings.simplefilter('ignore', chemicals.PlaceholderWarning)
    warnings.simplefilter('ignore', bst.exceptions.CostWarning)

    # 1. Reaction model: species, stoichiometry, rate law and parameters
    #    come from an EnzymeML document.
    document = enzymeml.read(examples.FDH_ENZYMEML)
    reaction = document.get_reaction('r_fdh')
    print(f"EnzymeML document: {document.name}")
    print(f"  rate law: {reaction.rate_law}")
    for parameter in document.parameters.values():
        print(f"  {parameter.symbol} = {parameter.value:g} "
              f"{parameter.unit_name}")

    # 2. Property package: chemicals are resolved from the identifiers of
    #    the document (InChIKey, ChEBI). The enzyme and NADH are not in the
    #    databases; their molecular weight and formula are supplied here.
    properties = {
        'nadh': {'formula': 'C21H29N7O14P2'},
        'fdh': {'MW': 40400.},
    }
    bst.main_flowsheet.set_flowsheet('fdh')
    bst.settings.set_thermo(document.create_chemicals(properties=properties))

    # 3. Flowsheet: feed pump, heater and stirred tank reactor.
    feed = bst.Stream('feed', Water=1000., units='kg/hr', T=298.15)
    feed.imol['formate'] = 0.1      # 100 mol/m3
    feed.imol['nad'] = 5e-3         # 5 mol/m3
    feed.imol['fdh'] = 2e-7         # 0.2 mmol/m3
    P1 = bst.Pump('P1', ins=feed, outs='pressurized_feed', P=2 * 101325)
    H1 = bst.HXutility(
        'H1', ins=P1-0, outs='heated_feed', T=document.temperature)
    R1 = bst.SinglePhaseReactor(
        'R1', ins=H1-0, outs='product',
        tau=2., T=document.temperature, P=101325,
    )
    document.parameterize(R1)  # Kinetics from EnzymeML
    fdh_sys = bst.System.from_units('fdh_sys', [P1, H1, R1])
    fdh_sys.simulate()
    product = R1.outs[0]
    print("\nSimulated in BioSTEAM:")
    print("  NAD+ conversion: "
          f"{1 - product.imol['nad'] / feed.imol['nad']:.2%}")
    print(f"  NADH produced:   {product.imol['nadh'] * 1e3:.3f} mol/hr")
    print(f"  reactor volume:  {R1.design_results['Reactor volume']:.2f} m3")

    # 4. Export the flowsheet as a DEXPI 2.0 process model.
    folder = folder or tempfile.mkdtemp()
    path = os.path.join(folder, 'fdh_cstr.xml')
    fdh_sys.to_dexpi(path)
    print("\nExported DEXPI process model "
          f"({os.path.getsize(path) / 1e3:.0f} kB)")

    # 5. Import it again into a separate flowsheet with a property package
    #    that is rebuilt from the identifiers stored in the document.
    copy = bst.System.from_dexpi(
        path, flowsheet='fdh_imported',
        resolver=chemicals.default_resolver(properties=properties),
    )
    copy.simulate()
    print("Re-imported from DEXPI:")
    print("  units:",
          ', '.join(f"{i.ID} ({type(i).__name__})" for i in copy.units))
    product_copy = copy.flowsheet.stream.product
    for ID in ('formate', 'nad', 'nadh', 'co2'):
        print(f"  {ID:8s} {product.imol[ID]:.9e} -> "
              f"{product_copy.imol[ID]:.9e} kmol/hr")
    difference = max(abs(product.mol - product_copy.mol))
    print(f"  largest difference in product flows: {difference:.1e} kmol/hr")
    return fdh_sys, copy


if __name__ == '__main__':
    main()
