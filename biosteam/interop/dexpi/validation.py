# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Structural validation of DEXPI process models.

The checks cover the rules of the DEXPI process information model that
BioSTEAM relies on to build a flowsheet: every stream connects an outlet
port to an inlet port, ports belong to process steps, identifiers are
present, and compositions match their material template.

"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from .model import ENGINEERING_MODEL, DexpiDocument, DexpiObject

__all__ = (
    'DexpiValidationError', 'ValidationIssue', 'find_issues', 'validate',
)


class DexpiValidationError(ValueError):
    """
    Raised when a document violates the DEXPI process information model.

    Attributes
    ----------
    issues : list[ValidationIssue]
        All issues that were found.

    """

    def __init__(self, issues):
        self.issues = list(issues)
        lines = '\n'.join(f' - {i.message}' for i in self.issues)
        super().__init__(f"invalid DEXPI process model:\n{lines}")


@dataclass(frozen=True)
class ValidationIssue:
    """A violation of the DEXPI process information model."""

    message: str
    obj: Optional[DexpiObject] = None


def _name(obj: DexpiObject) -> str:
    return f"{obj.class_name} '{obj.get('Identifier') or obj.id}'"


def _all_steps(steps):
    for step in steps:
        yield step
        yield from _all_steps(step.children('SubProcessSteps'))


def find_issues(document: DexpiDocument) -> list[ValidationIssue]:
    """
    Return all structural issues of a document.

    Examples
    --------
    >>> from biosteam.interop.dexpi import model, validation
    >>> root = model.DexpiObject(model.ENGINEERING_MODEL)
    >>> process = root.add(
    ...     'ConceptualModel', model.DexpiObject(model.PROCESS_MODEL))
    >>> stream = process.add('ProcessConnections', model.DexpiObject(
    ...     model.process_type('Stream'), data={'Identifier': 's1'}))
    >>> for issue in validation.find_issues(model.DexpiDocument(root)):
    ...     print(issue.message)
    Stream 's1' has no Source port
    Stream 's1' has no Target port

    """
    issues: list[ValidationIssue] = []

    def report(message, obj=None):
        issues.append(ValidationIssue(message, obj))

    if document.root.type != ENGINEERING_MODEL:
        report(f"root object is a '{document.root.type}', "
               "not an EngineeringModel",
               document.root)
    model = document.process_model
    if model is None:
        report("document does not contain a ProcessModel", document.root)
        return issues

    port_owner: dict[int, DexpiObject] = {}
    identifiers: dict[str, DexpiObject] = {}
    for step in _all_steps(model.children('ProcessSteps')):
        identifier = step.get('Identifier')
        if not identifier:
            report(f"{_name(step)} has no Identifier", step)
        elif identifier in identifiers:
            report(f"Identifier '{identifier}' is used by more than one "
                   "process step", step)
        else:
            identifiers[identifier] = step
        for port in step.children('Ports'):
            port_owner[id(port)] = step
            if not port.get('Identifier'):
                report(f"a port of {_name(step)} has no Identifier", port)
            if port.get('NominalDirection') is None:
                report(f"port '{port.get('Identifier')}' of {_name(step)} "
                       "has no NominalDirection", port)
            if port.target('ConnectorReference') is None:
                report(f"port '{port.get('Identifier')}' of {_name(step)} "
                       "is not connected", port)

    expected_direction = {'Source': 'Outlet', 'Target': 'Inlet'}
    for connection in model.children('ProcessConnections'):
        if not connection.get('Identifier'):
            report(f"{_name(connection)} has no Identifier", connection)
        for end, direction in expected_direction.items():
            port = connection.target(end)
            if port is None:
                report(f"{_name(connection)} has no {end} port", connection)
                continue
            if id(port) not in port_owner:
                report(f"the {end} of {_name(connection)} is not a port of "
                       "a process step", connection)
                continue
            actual = getattr(port.get('NominalDirection'), 'literal', None)
            if actual is not None and actual != direction:
                report(f"the {end} of {_name(connection)} is an {actual} "
                       f"port; expected an {direction} port", connection)
            reference = port.target('ConnectorReference')
            if reference is not None and reference is not connection:
                report(f"the {end} port of {_name(connection)} refers to "
                       f"{_name(reference)}", port)
        _check_composition(connection, report)
    return issues


def _check_composition(connection: DexpiObject, report) -> None:
    template = connection.target('MaterialTemplateReference')
    state = connection.target('MaterialStateReference')
    if template is None or state is None:
        return
    components = template.target('ListOfComponents')
    if components is None:
        report(f"MaterialTemplate '{template.get('Identifier')}' has no "
               "ListOfComponents", template)
        return
    size = len(components.targets('Component'))
    declared = template.get('NumberOfMaterialComponents')
    if isinstance(declared, int) and declared != size:
        report(f"MaterialTemplate '{template.get('Identifier')}' declares "
               f"{declared} components but lists {size}", template)
    for state_type in (*state.children('State'), *state.children('Phases')):
        composition = state_type.child('Composition')
        if composition is None:
            continue
        for name, values in composition.components.items():
            if len(values) not in (0, size):
                report(f"{name} of '{state_type.get('Identifier')}' has "
                       f"{len(values)} values for {size} components",
                       composition)


def validate(document: DexpiDocument) -> None:
    """
    Check the structure of a document.

    Raises
    ------
    DexpiValidationError
        If the document has structural issues.

    """
    issues = find_issues(document)
    if issues:
        raise DexpiValidationError(issues)
