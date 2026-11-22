# Supplementary data

Everything needed to reproduce the two use cases of the manuscript.

## Use case 1: NADH regeneration by formate dehydrogenase (stirred tank)

The input files ship with the software:

- `biosteam/interop/data/fdh_nadh_regeneration.json`: EnzymeML v2 document.
  The kinetic parameters are representative values, not measurements.
- `biosteam/interop/data/fdh_flowsheet.xml`: DEXPI 2.0 process model.

Run `docs/tutorial/Lab_data_to_simulation.ipynb`, or the tests in
`tests/test_fdh_cstr_roundtrip.py`, which assert every number reported.

## Use case 2: laccase flow reactor of Behr et al. (2024)

See `use_case_2_laccase/`.
