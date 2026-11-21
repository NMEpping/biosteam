# Design notes: DEXPI 2.0 and EnzymeML interoperability

This document records the design of `biosteam.interop`: what was audited,
which decisions were taken, and why. It is the reference for reviewers and
for anyone extending the mapping.

## 1. Scope

| Capability | Entry point |
|---|---|
| Export a `System` to a DEXPI 2.0 Process document | `System.to_dexpi(path)` / `interop.dexpi.to_dexpi(system)` |
| Rebuild an executable `System` from a DEXPI Process document | `System.from_dexpi(path)` / `interop.dexpi.from_dexpi(source)` |
| Resolve component identities to `thermosteam.Chemical` objects | `interop.chemicals` |
| Load EnzymeML reactions and rate laws as `thermosteam.Reaction` objects | `interop.enzymeml` |
| Reference case: FDH-catalysed NADH regeneration in a CSTR | `interop.examples`, `docs/tutorial/DEXPI_and_EnzymeML_interoperability.rst` |

Covered DEXPI content: the *Process* information model (`ProcessModel`,
`ProcessStep` and subclasses, `MaterialPort`, `Stream`, `MaterialTemplate`,
`ListOfMaterialComponents`, `Pure`/`CustomMaterialComponent`,
`MaterialState`, `Composition`, `QualifiedValue`, `Note`,
`PersistentIdentifier`). Out of scope for now: the *Plant* (P&ID) model,
diagram/graphics, instrumentation activities, and energy ports/flows (see §8).

## 2. Dependency audit

| Package | Version audited | License | Footprint | Python | Verdict |
|---|---|---|---|---|---|
| pyDEXPI | 2.1.0 (DEXPI 2.0.1 classes) | **AGPL-3.0** | lxml, networkx, openpyxl, pandas ≥ 3, pydantic ≥ 2.13 | ≥ 3.12 | Not a runtime dependency; optional *test-time* validator only |
| pyenzyme | 2.3.0 (EnzymeML v2) | BSD-2-Clause | bokeh, libsbml, mdmodels, rdflib, sympy, matplotlib, … | ≥ 3.11 | Optional (`biosteam[enzymeml]`), only needed for OMEX/SBML archives |

**pyDEXPI.** The Python versions line up with BioSTEAM (≥ 3.12) and the XML
serializer is faithful, but the license decides the question: BioSTEAM is
distributed under the UIUC/NCSA license, and code that imports an AGPL-3.0
library becomes subject to AGPL terms when the combination is conveyed or
offered as a network service. Making pyDEXPI a dependency — even an optional,
lazily imported one — would put that obligation on downstream users of a
permissively licensed package. pyDEXPI also pins recent `pandas`/`pydantic`
releases, which would constrain BioSTEAM environments.

The DEXPI 2.0 XML exchange format itself is small and regular (five element
kinds: `Object`, `Components`, `References`, `Data`, and typed values). The
decision is therefore:

* `biosteam.interop.dexpi` contains a **self-contained adapter written against
  the DEXPI specification** using only the standard library
  (`xml.etree.ElementTree`). No pyDEXPI code is imported or copied.
* pyDEXPI is used **only in the test suite**, as an independent conformance
  check: exported documents must load in pyDEXPI's strict mode, and documents
  written by pyDEXPI must import into BioSTEAM. Those tests are skipped when
  pyDEXPI is not installed (`pytest.importorskip`).

**pyenzyme.** The license is compatible, but the dependency tree is heavy for
what BioSTEAM needs (species, stoichiometry, a rate-law string and parameters).
EnzymeML v2 JSON is a plain JSON document, so it is parsed directly. pyenzyme
is imported lazily in exactly one place: reading `.omex` archives (which wrap
SBML and need libSBML). `pyenzyme.EnzymeMLDocument` objects are also accepted
by duck typing (`model_dump()`), without importing pyenzyme.

**Consequence:** `import biosteam` is unchanged — nothing in
`biosteam/__init__.py` imports `biosteam.interop`, and `biosteam.interop`
needs nothing beyond BioSTEAM's existing dependencies for DEXPI XML and
EnzymeML JSON. A missing optional package raises an `ImportError` naming the
extra to install, and only at the call that needs it.

## 3. Module layout

```
biosteam/interop/
├── __init__.py            lazy access to the submodules
├── chemicals.py           ChemicalIdentity + resolver chain (user data, thermosteam, PubChem, placeholder)
├── reactions.py           JSON codec for thermosteam reactions + registry of serializable rate laws
├── examples.py            FDH NADH-regeneration CSTR reference system
├── data/                  FDH EnzymeML dataset (JSON)
├── dexpi/
│   ├── model.py           DexpiObject tree and typed data values
│   ├── vocabulary.py      class hierarchy of DEXPI process steps
│   ├── xmlio.py           DEXPI XML reader / writer (stdlib ElementTree)
│   ├── quantities.py      DEXPI unit literals ⇄ BioSTEAM base units
│   ├── notes.py           BioSTEAM annotations as DEXPI notes
│   ├── values.py          JSON encoding of unit parameters
│   ├── registry.py        UnitCodec + UnitRegistry (BioSTEAM class ⇄ DEXPI class)
│   ├── codecs.py          codecs and default registrations for built-in units
│   ├── generic.py         GenericProcessStep fallback unit
│   ├── exporter.py        System → DexpiObject tree
│   ├── importer.py        DexpiObject tree → System
│   └── validation.py      structural validation of process models
└── enzymeml/
    ├── document.py        immutable records + JSON / OMEX / pyenzyme readers
    ├── units.py           EnzymeML UnitDefinition → (kmol, m, hr, kg, K)
    ├── ratelaw.py         safe rate-law evaluation (restricted `ast`)
    └── kinetics.py        SteadyStateCSTRRate, EnzymeML → thermosteam.Reaction
```

Each module has one reason to change: the wire format (`xmlio`), the vocabulary
(`model`, `vocabulary`), unit conventions (`quantities`), class mapping (`registry`,
`codecs`), and the two directions of translation (`exporter`, `importer`).

## 4. DEXPI mapping

### 4.1 Structure

| BioSTEAM | DEXPI 2.0 Process |
|---|---|
| `System` | `EngineeringModel` → `ProcessModel` |
| `Unit` | `ProcessStep` subclass chosen by the registry |
| inlet / outlet slot | `MaterialPort` (`NominalDirection` Inlet/Outlet), identifier `in<i>` / `out<i>` |
| `Stream` | `Stream` with `Source`/`Target` port references |
| feed / product stream | additional `Source` / `Sink` process step (DEXPI requires both stream ends) |
| `Chemicals` | `MaterialTemplate` + `ListOfMaterialComponents` |
| `Chemical` found in a database | `PureMaterialComponent` (+ `PersistentIdentifier`s for CAS, InChIKey, PubChem) |
| user-defined `Chemical` | `CustomMaterialComponent` |
| `Stream.T`, `.P`, `.F_mass` | `Stream.Temperature` [K], `.Pressure` [Pa], `.MassFlow` [kg/h] |
| `Stream.imol` | `MaterialState` → per-phase `MaterialStateType.Composition.MoleFlow` [kmol/s] |

### 4.2 Unit registry

`UnitRegistry` maps BioSTEAM classes to DEXPI classes through `UnitCodec`
strategy objects. Lookup walks the MRO on export (a subclass of `Pump` is a
`Pumping` step unless registered otherwise) and the DEXPI class hierarchy on
import (`RotaryMixing` → `Mixing` → `Mixer`). Default registrations:

| BioSTEAM | DEXPI | Import default for |
|---|---|---|
| `SinglePhaseReactor` (and other `AbstractStirredTankReactor`s) | `ReactingChemicals` (`Method=Tank`) | `ReactingChemicals` |
| `Mixer` | `MixingSimple` | `Mixing` |
| `Splitter` | `SplittingMaterial` | `SplittingMaterial` |
| `Pump` | `Pumping` | `Pumping` |
| `HXutility` | `ExchangingThermalEnergy` | `ExchangingThermalEnergy`, `SupplyingThermalEnergy`, `RemovingThermalEnergy` |
| `Flash` | `SeparatingByFlash` | `SeparatingByFlash` |
| `StorageTank` / `MixTank` | `StoringInTank` | `StoringInTank` |
| `IsenthalpicValve` | `RegulatingFlow` | `RegulatingFlow` |
| any other `Unit` | `ProcessStep` | — |
| — | any unmapped class | `GenericProcessStep` |

Third parties extend the mapping without touching BioSTEAM:

```python
from biosteam.interop import dexpi
dexpi.default_registry.register(MyCentrifuge, 'SeparatingByCentrifugalForce',
                                import_default=True)
```

### 4.3 BioSTEAM-specific data: `Note` annotations

DEXPI has no slot for simulator-specific parameters (splits, reactions,
residence times, costing options). Rather than inventing non-standard XML,
every exported object that needs them references a standard `Core/Note` with
`NoteClassification = "BioSTEAM"` whose text is a JSON payload:

```json
{"schema": 1, "kind": "unit", "class": "SinglePhaseReactor",
 "parameters": {"tau": 2.0, "T": 303.15, "reactions": {...}}}
```

Any DEXPI tool can read, display and re-export notes, so the payload survives
third-party round trips, and a document remains valid DEXPI for tools that
ignore it. Parameters are collected generically from the unit's `_init`
signature (BioSTEAM's convention is that constructor arguments are stored as
same-named attributes); codecs override this where a parameter is not plain
data (splits, reactions).

### 4.4 Never losing topology

* **Unmapped BioSTEAM unit → DEXPI:** exported as a plain `ProcessStep` with
  its ports, streams and a note carrying class name and parameters.
* **Unmapped DEXPI step → BioSTEAM:** imported as `GenericProcessStep`, a
  pass-through unit that remembers the original DEXPI class and data
  attributes and writes them back on export.
* **Annotated class not available on import:** falls back to the registry
  default for the DEXPI class, then to `GenericProcessStep`, with a warning.

Class names from a document are resolved **only** against registered classes
and the unit classes in the `biosteam` namespace (BioSTEAM adds user-defined
`Unit` subclasses there when they are created, so custom units round-trip
within a session that defines them) — never by importing a dotted path — so
a document cannot make BioSTEAM import arbitrary modules. The XML reader
refuses documents with DTD or entity declarations.

## 5. Chemical identity and property hydration

Topology and thermodynamics are deliberately separate. A document yields
`ChemicalIdentity` records (ID, name, CAS, InChIKey, InChI, SMILES, ChEBI,
PubChem CID, formula, MW, phase); a `ResolverChain` turns them into
`thermosteam.Chemical` objects:

1. `UserChemicalsResolver` — chemicals the caller already has (matched by ID,
   CAS or InChIKey).
2. `PropertyDatabaseResolver` — user-supplied property records.
3. `ThermosteamResolver` — thermosteam's databases, tried by CAS, InChIKey,
   PubChem CID, SMILES, name, then ID.
4. `PubChemResolver` — optional, off by default; uses PUG-REST to obtain
   formula/MW/InChIKey (ChEBI IDs are resolved through PubChem synonyms).
5. `PlaceholderResolver` — last resort: a database-free chemical from the
   known formula/MW with default properties, emitted with a warning.

Chemicals that a database knows only partially (NAD⁺ has a formula but no
property models) are fixed to one phase and completed with water-like
defaults, again with a warning. The molar heat capacity of macromolecules is
capped, because water's specific heat times an enzyme's molecular weight
exceeds the range the property models accept. thermosteam has no field for
ChEBI identifiers; they are kept as chemical aliases (`'CHEBI:15846'`).

Both the DEXPI importer and the EnzymeML adapter produce `ChemicalIdentity`
objects, so species from either source align on the same identifiers.
`from_dexpi(..., chemicals=...)` or `resolver=...` lets the caller inject a
property package; if the active thermo already defines every component it is
used as is.

## 6. EnzymeML bridge

`enzymeml.read(source)` accepts a `.json` path, a `.omex` path (pyenzyme),
a `dict`, or a pyenzyme document, and returns an immutable
`EnzymeMLDocument`. From it:

* `document.chemical_identities()` → identities for §5,
* `document.create_chemicals(resolver)` → `thermosteam.Chemicals`,
* `document.to_reaction(reaction_id, residence_time=...)` →
  `thermosteam.Reaction` with a kinetic `rate`.

**Rate laws.** The EnzymeML equation string is compiled by a restricted
`ast` walker (arithmetic, powers, and a whitelist of math functions; no
attribute access, calls to anything else, or names outside the species and
parameter tables). Parameters are converted with their EnzymeML
`UnitDefinition` to a single coherent system (kmol, m³, h), so evaluating the
law at concentrations in kmol/m³ yields kmol/(m³·h).

**CSTR semantics.** thermosteam evaluates `Reaction.rate(stream)` once and
adds `rate × stoichiometry` to the stream. `SteadyStateCSTRRate` implements
that callable as the solution of the CSTR design equation
`ξ = r(C_out(ξ)) · V`, with `C_out,i = (F_i + ν_i ξ) / Q` and `V = τ · Q`,
solved by bracketed root finding (Brent's method). BioSTEAM has no unit
named `CSTR`; its continuous stirred tank is `SinglePhaseReactor`, and the
rate works inside it without a new unit class, reads `tau` live from the
reactor it is bound to, and is serializable (`to_dict` / `from_dict`), which
is how kinetics travel inside the DEXPI note and are restored on import.

## 7. Testing strategy

* `tests/test_dexpi.py` — XML reader/writer, structural validation, registry
  lookup, generic fallbacks, multi-unit and recycle round trips, pyDEXPI
  conformance (skipped without pyDEXPI).
* `tests/test_enzymeml.py` — parsing, unit conversion, rate-law safety,
  CSTR balance against the analytical solution, missing-dependency handling.
* `tests/test_interop_chemicals.py` — resolver chain and fallbacks.
* `tests/test_fdh_cstr_roundtrip.py` — EnzymeML → simulate → DEXPI →
  re-import in a separate flowsheet with a rebuilt property package →
  identical stream, design and cost results; a `slow` test repeats this
  across two Python processes that share only the DEXPI file.

## 8. Known limitations and next steps

* Energy ports/flows (`ThermalEnergyPort`, `ThermalEnergyFlow`) are not
  emitted; utility duties stay in BioSTEAM results.
* Nested `ProcessStep`s (`SubProcessSteps`) are flattened on import;
  subsystems are exported flat.
* Process specifications, custom `_run` overrides and facilities are not
  serializable; affected units are exported with a warning in the note.
* DEXPI's mole-flow unit list has no kmol/h, so compositions are written in
  kmol/s; stored stream values agree to floating-point precision (~1e-15
  relative), not necessarily bitwise.
* A re-imported system reproduces the original results when simulated from
  the same starting point. Systems with recycles converge within the
  system's tolerance, so results can differ at that level if the imported
  stream values are used as the initial guess.
* Only the DEXPI XML format is read and written. pyDEXPI's JSON format is a
  format of that library, not of the DEXPI specification.
* Enzymes are modelled as chemicals on a molar basis with water-like
  properties; their energy-balance contribution is negligible at catalytic
  concentrations but is an estimate.
* Rate laws are isothermal parameter sets; temperature-dependent parameters
  would need an Arrhenius extension of the EnzymeML document.
