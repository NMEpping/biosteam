# Figures

Each figure is a standalone TikZ document that compiles to a vector PDF:

    cd figures && for f in fig_*.tex; do pdflatex $f; done

| File | Content |
|---|---|
| `fig_workflow.pdf` | Data flow and FAIR aspects (also suitable as graphical abstract) |
| `fig_routes_levels.pdf` | (a) Ontology-based and format-based route across the levels ontology, data model and exchange format; (b) machine-readable, -interpretable, -actionable |
| `fig_laccase.pdf` | Use case 2; written by `supplementary/use_case_2_laccase/run.py` (matplotlib, same fonts and colours) |

Design (`tuhhstyle.tex`): TUHH colour palette; Latin Modern as in the
manuscript; highlights bold at caption size (10 pt), longer text regular
at 8 pt. The figures are drawn at the manuscript's text width (390 pt)
and must be included unscaled so that these sizes hold.

## Logos

`\logo{name}{height}{fallback}` uses `logos/<name>.pdf` or `.png` if the
file exists and typesets the name as a wordmark otherwise.

- `logos/biosteam.png`: from this repository (`docs/images/logo.png`).
- `logos/python.pdf`: the Python logo as distributed with `ipykernel`.
  It is a trademark of the Python Software Foundation; see
  https://www.python.org/psf/trademarks/ for the usage policy.
- `logos/enzymeml.pdf` and `logos/dexpi.pdf` are **not included**: the
  official files could not be obtained when the figures were made, so both
  currently appear as wordmarks. Add the official logo files under these
  names and recompile. Check with the EnzymeML team and DEXPI e.V. that
  use in a publication is permitted.
