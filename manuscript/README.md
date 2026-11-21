# Short communication (draft)

*Formats carry, ontologies anchor: FAIR exchange of enzyme kinetics and
process flowsheets between EnzymeML, DEXPI and BioSTEAM.*

- Target journal: *Digital Chemical Engineering* (Elsevier), as a short
  communication. Check the journal's guide for authors for the current
  article types and length limits before submission.
- Template: Elsevier `elsarticle` class with the `elsarticle-num`
  reference style (both part of TeX Live, package `texlive-publishers`).

Build:

    pdflatex manuscript && bibtex manuscript && pdflatex manuscript && pdflatex manuscript

Before submission, complete the items marked `TODO (authors)` in
`manuscript.tex`: author names and affiliations, CRediT statement, funding,
and an archived release with a DOI for the data availability statement.

This directory is not part of the BioSTEAM package and should be left out
of a pull request to the upstream repository.
