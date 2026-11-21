DEXPI and EnzymeML interoperability
===================================

BioSTEAM can exchange process flowsheets with other engineering tools through
`DEXPI <https://dexpi.org>`_ (Data Exchange in the Process Industry) and load
enzyme kinetics from `EnzymeML <https://enzymeml.org>`_ documents. This
tutorial follows one case through the complete chain:

1. load a reaction model from EnzymeML,
2. resolve its species to chemicals,
3. simulate a continuous stirred tank reactor (CSTR),
4. export the flowsheet as a DEXPI 2.0 process model, and
5. import it again and reproduce the results.

The case is the regeneration of the cofactor NADH by formate dehydrogenase
(FDH), with formate as a sacrificial cosubstrate:

.. math:: \text{Formate} + \text{NAD}^+ \xrightarrow{\text{FDH}} \text{CO}_2 + \text{NADH}

The complete script is :download:`dexpi_fdh_cstr.py`. Neither DEXPI XML nor
EnzymeML JSON needs packages beyond BioSTEAM itself.

Reaction model from EnzymeML
----------------------------

An EnzymeML document records the species of a reaction, its stoichiometry,
a rate law and kinetic parameters. BioSTEAM ships the document of the
reference case:

.. code-block:: python

    >>> import biosteam as bst
    >>> from biosteam.interop import enzymeml, examples
    >>> document = enzymeml.read(examples.FDH_ENZYMEML)
    >>> reaction = document.get_reaction('r_fdh')
    >>> reaction.rate_law
    'kcat * fdh * nad * formate / (K_i_nad * K_m_formate + K_m_formate * nad + K_m_nad * formate + nad * formate)'
    >>> for parameter in document.parameters.values():
    ...     print(parameter.symbol, parameter.value, parameter.unit_name)
    kcat 3.7 1 / s
    K_m_formate 5.6 mmol / l
    K_m_nad 0.045 mmol / l
    K_i_nad 0.04 mmol / l

The rate law is the initial-rate form of an ordered bi-bi mechanism. The
parameters of this document are representative values for illustration, not
a curated measurement; replace the file with your own data for real work.
:func:`~biosteam.interop.enzymeml.read` also accepts dictionaries,
``pyenzyme`` documents and, with ``pip install biosteam[enzymeml]``, OMEX
archives.

Chemicals from identifiers
--------------------------

The document identifies its species by InChIKey and ChEBI number. These are
resolved to chemicals in a fixed order: chemicals you provide, property
records you provide, thermosteam's databases, optionally PubChem
(``online=True``), and finally placeholders with water-like properties.
Formate, NAD\ :sup:`+` and CO\ :sub:`2` are found in the databases; the
formula of NADH and the molecular weight of the enzyme are supplied as
property records:

.. code-block:: python

    >>> properties = {
    ...     'nadh': {'formula': 'C21H29N7O14P2'},
    ...     'fdh': {'MW': 40400.},
    ... }
    >>> chemicals = document.create_chemicals(properties=properties)
    >>> chemicals.IDs
    ('Water', 'formate', 'nad', 'co2', 'nadh', 'fdh')
    >>> chemicals.formate.CAS
    '64-18-6'
    >>> bst.settings.set_thermo(chemicals)

.. note::

    Properties that no source provides (for example the heat capacity of
    NAD\ :sup:`+`) default to those of water, and a ``PlaceholderWarning``
    names them. Mass balances are unaffected; treat energy balances of such
    species as estimates.

Simulating the CSTR
-------------------

The reactor is an ordinary :class:`~biosteam.units.SinglePhaseReactor`.
:meth:`~biosteam.interop.enzymeml.EnzymeMLDocument.parameterize` assigns the
reaction and ties its kinetics to the residence time of the reactor:

.. code-block:: python

    >>> feed = bst.Stream('feed', Water=1000., units='kg/hr', T=298.15)
    >>> feed.imol['formate'] = 0.1   # 100 mol/m3
    >>> feed.imol['nad'] = 5e-3      # 5 mol/m3
    >>> feed.imol['fdh'] = 2e-7      # 0.2 mmol/m3
    >>> P1 = bst.Pump('P1', ins=feed, outs='pressurized_feed', P=2 * 101325)
    >>> H1 = bst.HXutility('H1', ins=P1-0, outs='heated_feed', T=document.temperature)
    >>> R1 = bst.SinglePhaseReactor('R1', ins=H1-0, outs='product',
    ...                             tau=2., T=document.temperature, P=101325)
    >>> reaction = document.parameterize(R1)
    >>> fdh_sys = bst.System.from_units('fdh_sys', [P1, H1, R1])
    >>> fdh_sys.simulate()
    >>> product = R1.outs[0]
    >>> print(f"{1 - product.imol['nad'] / feed.imol['nad']:.2%}")
    91.21%

The conversion follows from the steady-state balance of an ideal CSTR, in
which the reaction proceeds at outlet concentrations:
:math:`\xi = r(C_{out}) \cdot V` with :math:`V = \tau Q`. Changing ``R1.tau``
or the enzyme concentration in the feed changes the conversion accordingly.

Export to DEXPI
---------------

.. code-block:: python

    >>> fdh_sys.to_dexpi('fdh_cstr.xml')

The file is a DEXPI 2.0 process model. Each unit operation is a process
step (``Pumping``, ``ExchangingThermalEnergy``, ``ReactingChemicals``), each
stream a DEXPI ``Stream`` with temperature, pressure, flow and composition,
and each chemical a material component with its CAS number, InChIKey and
ChEBI identifier. Feeds and products end in ``Source`` and ``Sink`` steps.

DEXPI has no attributes for simulator-specific parameters such as residence
times or kinetics. BioSTEAM stores them as JSON in standard DEXPI ``Note``
objects, so the file stays valid for other tools and nothing is lost.

Import from DEXPI
-----------------

.. code-block:: python

    >>> from biosteam.interop import chemicals as chems
    >>> copy = bst.System.from_dexpi(
    ...     'fdh_cstr.xml', flowsheet='fdh_imported',
    ...     resolver=chems.default_resolver(properties=properties),
    ... )
    >>> copy.simulate()
    >>> [type(i).__name__ for i in copy.units]
    ['Pump', 'HXutility', 'SinglePhaseReactor']
    >>> product_copy = copy.flowsheet.stream.product
    >>> abs(product_copy.imol['nadh'] - product.imol['nadh']) < 1e-12
    True

The imported system lives in its own flowsheet and uses a property package
rebuilt from the identifiers in the file. Pass ``chemicals=`` instead of a
resolver to use a property package of your own.

Documents from other tools
--------------------------

Process steps without BioSTEAM annotations are created from the unit class
registered for their DEXPI class, using the standard attributes they carry
(for example the pressure of a ``Pumping`` step). Steps that BioSTEAM has no
model for become :class:`~biosteam.interop.dexpi.GenericProcessStep` units,
which keep the flowsheet connected and write their original class and
attributes back on export. Register your own classes to extend the mapping:

.. code-block:: python

    >>> from biosteam.interop import dexpi
    >>> class Centrifuge(bst.Splitter): pass
    >>> dexpi.default_registry.register(
    ...     Centrifuge, 'SeparatingByCentrifugalForce', import_default=True)

See :doc:`../API/interop/index` for the reference and the design notes.
