"""This module defines the Interpolate class."""

# Copyright (C) 2021 Nacime Bouziani
#
# This file is part of UFL (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Modified by Nacime Bouziani, 2021-2022

from ufl.argument import Argument, Coargument
from ufl.constantvalue import as_ufl
from ufl.core.base_form_operator import BaseFormOperator
from ufl.core.ufl_type import ufl_type
from ufl.duals import is_dual
from ufl.form import BaseForm
from ufl.functionspace import AbstractFunctionSpace


@ufl_type(num_ops="varying", is_differential=True)
class Interpolate(BaseFormOperator):
    """Symbolic representation of the interpolation operator."""

    # Slots are disabled here because they cause trouble in PyDOLFIN
    # multiple inheritance pattern:
    _ufl_noslots_ = True

    def __init__(self, expr, v):
        """Initialise.

        Args:
            expr: a UFL expression to interpolate.
            v: the FunctionSpace to interpolate into or the Coargument
                defined on the dual of the FunctionSpace to interpolate into.
        """
        dual_args = (Coargument, BaseForm)

        expr = as_ufl(expr)
        if isinstance(expr, dual_args):
            raise ValueError("Expecting the first argument to be primal.")

        if isinstance(v, AbstractFunctionSpace):
            if is_dual(v):
                raise ValueError("Expecting a primal function space.")
            from ufl.algorithms import extract_arguments

            expr_args = extract_arguments(expr)
            is_adjoint = len(expr_args) and expr_args[0].number() == 0
            v = Argument(v.dual(), 1 if is_adjoint else 0)
        elif not isinstance(v, dual_args):
            raise ValueError(
                "Expecting the second argument to be FunctionSpace, Coargument, or BaseForm."
            )

        # Reversed order convention
        argument_slots = (v, expr)
        # Get the primal space (V** = V)
        arg, *_ = v.arguments()
        function_space = arg.ufl_function_space()

        # Set the operand as `expr` for DAG traversal purpose.
        operand = expr
        BaseFormOperator.__init__(
            self, operand, function_space=function_space, argument_slots=argument_slots
        )

    def _ufl_expr_reconstruct_(self, expr, v=None, **add_kwargs):
        """Return a new object of the same type with new operands."""
        v = v or self.argument_slots()[0]
        return type(self)(expr, v, **add_kwargs)

    def __repr__(self):
        """Default repr string construction for Interpolate."""
        r = "Interpolate("
        r += ", ".join(repr(arg) for arg in reversed(self.argument_slots()))
        r += f"; {self.ufl_function_space()!r})"
        return r

    def __str__(self):
        """Default str string construction for Interpolate."""
        s = "Interpolate("
        s += ", ".join(str(arg) for arg in reversed(self.argument_slots()))
        s += f"; {self.ufl_function_space()})"
        return s

    def __eq__(self, other):
        """Check for equality."""
        if self is other:
            return True
        return (
            type(self) is type(other)
            and all(a == b for a, b in zip(self._argument_slots, other._argument_slots))
            and self.ufl_function_space() == other.ufl_function_space()
        )


# Helper function
def interpolate(expr, v):
    """Create symbolic representation of the interpolation operator.

    Args:
        expr: a UFL expression to interpolate.
        v: the FunctionSpace to interpolate into or the Coargument
            defined on the dual of the FunctionSpace to interpolate into.
    """
    return Interpolate(expr, v)
