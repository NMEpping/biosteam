# Use case 2: laccase-catalysed oxidation of ABTS in a capillary flow reactor

Reuse of the data published with

> Behr, A.S.; Surkamp, J.; Abbaspour, E.; Häußler, M.; Lütz, S.; Pleiss, J.;
> Kockmann, N.; Rosenthal, K. Fluent Integration of Laboratory Data into
> Biocatalytic Process Simulation Using EnzymeML, DWSIM, and Ontologies.
> *Processes* 2024, 12, 597. https://doi.org/10.3390/pr12030597

Source of the data: https://github.com/TUDoAD/DWSIM-EnzymeML-KG
(commit `2787ea7`, 2024-10-25), https://doi.org/10.5281/zenodo.10613473,
licensed under **CC BY 4.0**. The files `laccase_abts.enzymeml.json` and
`laccase_abts.dexpi.xml` are adaptations of that data (format conversion)
and are shared under the same licence.

| File | Content |
|---|---|
| `build_inputs.py` | Converts the two spreadsheets (`ELNs/`) of the source repository |
| `laccase_abts.enzymeml.json` | The EnzymeML spreadsheet as EnzymeML v2 JSON (2 vessels, 3 reactions, 12 measurements) |
| `laccase_abts.dexpi.xml` | The process flow diagram of the process spreadsheet as a DEXPI 2.0 process model |
| `run.py` | Simulation in BioSTEAM, comparison with the measurements, DEXPI round trip |
| `results.csv`, `fig_laccase.pdf`, `laccase_abts_simulated.dexpi.xml` | Output of `run.py` |
| `laccase_abts_diagram.*` | Flowsheet with a DEXPI diagram and its drawing (`../render_dexpi.py`) |

Reproduce:

    git clone https://github.com/TUDoAD/DWSIM-EnzymeML-KG
    python build_inputs.py DWSIM-EnzymeML-KG
    python run.py

Both generated files were checked with independent implementations: the
EnzymeML document loads with `pyenzyme` 2.3.0 and the DEXPI file in the
strict mode of `pyDEXPI` 2.1.0.

## Choices made in the conversion

No values were added or changed. Where the target formats ask for
something the spreadsheets state differently:

1. **Stoichiometry.** The EnzymeML spreadsheet lists educts and products
   without numbers; the process spreadsheet gives
   4 ABTS_red + O2 -> 4 ABTS_ox + 2 H2O. The EnzymeML document uses these
   numbers divided by four (per mole of ABTS), because the rate law
   `kcat * E * S / (Km + S)` gives the rate of ABTS oxidation.
2. **Water** is added as species `s6`, since it takes part in the reaction
   according to the process spreadsheet but is not a species of the
   EnzymeML spreadsheet.
3. **Species per vessel.** The spreadsheet lists each compound once per
   vessel (`s0`-`s2` straight, `s3`-`s5` helical capillary). Reactions `r0`
   and `r1` use the species of their vessel. Rate-law equations are copied
   verbatim, so the equation of `r0` refers to `s0` as in the original.
4. **Temperature** is converted from 38 °C to 311.15 K in the EnzymeML
   document and kept in °C in the DEXPI file.
5. **Enzyme concentration** of each experiment is stored as the initial
   value of the protein species.
6. **Flowsheet.** Material streams, the mixer and the plug-flow reactor of
   the `PFD` sheet become DEXPI streams, a `MixingSimple` and a
   `ReactingChemicals` step with `Method = Tubular`. Source and sink steps
   terminate feeds and product, as DEXPI requires. The energy stream
   `Heat` is not represented. The pressure drop, length and diameter of the
   reactor have no DEXPI attribute and are not transferred; volume,
   pressure and temperature are.
7. **CAS numbers** are taken from the process spreadsheet for ABTS and
   added for water and oxygen, which the spreadsheet marks as database
   compounds.

## What `run.py` needs beyond the two files

- a table assigning the EnzymeML species to the components of the
  flowsheet (`CHEMICAL_IDS`);
- the molecular weights of ABTS (reduced, oxidised) and laccase, as stated
  in the process spreadsheet (`PROPERTIES`);
- one derivation: residence time = reactor volume / sum of the volumetric
  feed flows, both read from the DEXPI file (60 s).

## Results

- Plug-flow simulation in BioSTEAM: 0.495 mM oxidised ABTS at 60 s,
  matching the integrated Michaelis-Menten equation to within 1e-5 mM.
- DWSIM result of Behr et al., from their knowledge graph
  (`ontologies/KG-DWSIM_EnzML_ELN_output.owl`): 0.448 mM. The 10 %
  difference is not explained; DWSIM was not rerun here.
- Measured (11.6 µM laccase, 0.97 mM ABTS, 1.34 mM oxygen): 0.32 mM at
  61 s in the helically coiled and 0.28 mM at 63 s in the straight
  capillary. Both simulations overpredict, as discussed by Behr et al.
- Only mass and energy balances are computed. BioSTEAM's vessel design
  and costing are made for stirred tanks and fail for a capillary of
  8 mL. The inlet concentrations are 0.962 mM ABTS and 11.5 µM laccase,
  because BioSTEAM computes the liquid volume from the molar flows at
  38 °C rather than taking the stated 5 mL/min.
