# Supplementary data

Everything needed to reproduce the two use cases of the manuscript. Each
folder holds the EnzymeML document, the DEXPI process model, a script that
runs the case, the simulated system exported as DEXPI, and a drawing of
the flowsheet.

| | Use case 1 (`use_case_1_fdh/`) | Use case 2 (`use_case_2_laccase/`) |
|---|---|---|
| System | NADH regeneration by formate dehydrogenase, stirred tank | Laccase-catalysed oxidation of ABTS, capillary flow reactor |
| Data | Illustrative parameters (not measurements) | Behr et al. (2024), CC BY 4.0, converted by `build_inputs.py` |
| EnzymeML v2 JSON | `fdh_nadh_regeneration.enzymeml.json` | `laccase_abts.enzymeml.json` |
| DEXPI 2.0 flowsheet | `fdh_flowsheet.dexpi.xml` | `laccase_abts.dexpi.xml` |
| Run | `run.py` | `run.py` |
| Simulated system as DEXPI | `fdh_simulated.dexpi.xml` | `laccase_abts_simulated.dexpi.xml` |
| Flowsheet with DEXPI diagram | `fdh_flowsheet_diagram.dexpi.xml` | `laccase_abts_diagram.dexpi.xml` |
| Drawing | `fdh_flowsheet_diagram.svg`, `.pdf` | `laccase_abts_diagram.svg`, `.pdf` |

The two input files of use case 1 are copies of
`biosteam/interop/data/fdh_nadh_regeneration.json` and `fdh_flowsheet.xml`.

## Drawings (`render_dexpi.py`)

The flowsheets are drawn with pyDEXPI 2.1.0:

    python render_dexpi.py use_case_1_fdh/fdh_flowsheet.dexpi.xml
    python render_dexpi.py use_case_2_laccase/laccase_abts.dexpi.xml

A DEXPI process model does not have to contain a drawing, and pyDEXPI
renders drawings but does not lay them out. The script therefore places
the process steps from left to right, adds the result to the model as a
DEXPI `Diagram` (blocks, port positions, connector lines, flow arrows,
labels) using pyDEXPI's classes, saves that model with pyDEXPI's XML
serializer, and renders it with pyDEXPI's SVG exporter. The PDF is
converted from the SVG with `cairosvg`. Blocks show the identifier and
the DEXPI class of each step.

The files with a diagram still load as process models: in pyDEXPI's
strict mode and in BioSTEAM, which ignores the drawing.

pyDEXPI is licensed under the AGPL-3.0. `render_dexpi.py` imports it and
is a stand-alone tool; it is not part of BioSTEAM, which does not depend
on pyDEXPI.
