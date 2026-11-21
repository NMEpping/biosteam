# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2026, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Safe evaluation of rate law expressions.

Rate laws arrive as text from documents that may not be trusted. They are
parsed into a syntax tree and evaluated by walking the tree; only numbers,
names, arithmetic and a fixed set of mathematical functions are accepted,
so that an expression can never execute arbitrary code.

"""
from __future__ import annotations

import ast
import math
import operator
from typing import Callable, Iterable, Mapping, Optional

__all__ = ('RateLaw', 'RateLawError')

_BINARY_OPERATORS = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.Pow: operator.pow,
}

_UNARY_OPERATORS = {
    ast.USub: operator.neg,
    ast.UAdd: operator.pos,
}

#: Functions that may be called in a rate law.
FUNCTIONS: dict[str, Callable] = {
    'exp': math.exp,
    'log': math.log,
    'ln': math.log,
    'log10': math.log10,
    'sqrt': math.sqrt,
    'abs': abs,
    'min': min,
    'max': max,
    'pow': math.pow,
    'tanh': math.tanh,
}

#: Longest expression accepted [characters].
MAX_LENGTH = 10_000


class RateLawError(ValueError):
    """Raised when a rate law is malformed or cannot be evaluated."""


class RateLaw:
    """
    A rate law expression that can be evaluated for given variable values.

    Parameters
    ----------
    equation :
        Expression in infix notation. Both ``**`` and ``^`` denote powers.
    names :
        Names that may appear in the expression. If given, any other name
        is rejected.

    Examples
    --------
    >>> from biosteam.interop.enzymeml.ratelaw import RateLaw
    >>> law = RateLaw('kcat * E * S / (K_m + S)')
    >>> sorted(law.symbols)
    ['E', 'K_m', 'S', 'kcat']
    >>> law({'kcat': 10., 'E': 0.5, 'S': 2., 'K_m': 2.})
    2.5

    Anything but arithmetic is rejected:

    >>> RateLaw('__import__("os").system("ls")')
    Traceback (most recent call last):
    RateLawError: function '__import__' is not allowed in a rate law

    """

    __slots__ = ('equation', 'symbols', '_tree')

    def __init__(self, equation: str, names: Optional[Iterable[str]] = None):
        if not isinstance(equation, str) or not equation.strip():
            raise RateLawError("rate law must be a non-empty string")
        if len(equation) > MAX_LENGTH:
            raise RateLawError("rate law is too long")
        self.equation = equation
        try:
            tree = ast.parse(equation.strip().replace('^', '**'), mode='eval')
        except SyntaxError as error:
            raise RateLawError(
                f"rate law {equation!r} is not a valid expression: {error.msg}"
            ) from None
        self._tree = tree.body
        self.symbols = frozenset(self._validate(self._tree))
        if names is not None:
            unknown = sorted(self.symbols.difference(names))
            if unknown:
                raise RateLawError(
                    f"rate law {equation!r} uses undefined names {unknown}"
                )

    def _validate(self, node) -> set[str]:
        """Return the names used below `node`, rejecting unsafe syntax."""
        if isinstance(node, ast.Constant):
            if (isinstance(node.value, bool)
                    or not isinstance(node.value, (int, float))):
                raise RateLawError(
                    "only numbers are allowed as constants in a rate law")
            return set()
        if isinstance(node, ast.Name):
            return {node.id}
        if isinstance(node, ast.BinOp) and type(node.op) in _BINARY_OPERATORS:
            return self._validate(node.left) | self._validate(node.right)
        if isinstance(node, ast.UnaryOp) and type(node.op) in _UNARY_OPERATORS:
            return self._validate(node.operand)
        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name) or node.keywords:
                raise RateLawError(
                    "only plain function calls are allowed in a rate law")
            if node.func.id not in FUNCTIONS:
                raise RateLawError(
                    f"function '{node.func.id}' is not allowed in a rate law"
                )
            names = set()
            for argument in node.args:
                names |= self._validate(argument)
            return names
        raise RateLawError(
            f"'{type(node).__name__}' syntax is not allowed in a rate law"
        )

    def _evaluate(self, node, values: Mapping[str, float]) -> float:
        if isinstance(node, ast.Constant):
            # Floats keep powers of large constants from growing unboundedly.
            return float(node.value)
        if isinstance(node, ast.Name):
            return values[node.id]
        if isinstance(node, ast.BinOp):
            return _BINARY_OPERATORS[type(node.op)](
                self._evaluate(node.left, values),
                self._evaluate(node.right, values),
            )
        if isinstance(node, ast.UnaryOp):
            return _UNARY_OPERATORS[type(node.op)](
                self._evaluate(node.operand, values)
            )
        arguments = [self._evaluate(i, values) for i in node.args]
        return FUNCTIONS[node.func.id](*arguments)

    def __call__(self, values: Mapping[str, float]) -> float:
        """
        Return the value of the rate law.

        Parameters
        ----------
        values :
            Value of each name used in the expression.

        Raises
        ------
        RateLawError
            If a name has no value or the expression cannot be computed.

        """
        try:
            return float(self._evaluate(self._tree, values))
        except KeyError as error:
            raise RateLawError(
                f"no value for '{error.args[0]}' in rate law {self.equation!r}"
            ) from None
        except (ArithmeticError, ValueError, TypeError) as error:
            raise RateLawError(
                f"rate law {self.equation!r} could not be evaluated: {error}"
            ) from None

    def __repr__(self):
        return f"RateLaw({self.equation!r})"
