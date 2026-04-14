"""Algorithm for replacing gradients in an expression."""

# Copyright (C) 2008-2016 Martin Sandve Alnæs
#
# This file is part of UFL (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import functools

import ufl.classes
from ufl.algorithms.map_integrands import map_integrand_dags
from ufl.classes import ReferenceValue
from ufl.corealg.dag_traverser import DAGTraverser


class FunctionPullbackApplier(DAGTraverser):
    """A pull back applier."""

    @functools.singledispatchmethod
    def process(self, o: ufl.classes.Expr) -> ufl.classes.Expr:
        """Process ``o``."""
        return super().process(o)

    @process.register(ufl.classes.Expr)
    def _(self, o: ufl.classes.Expr) -> ufl.classes.Expr:
        return self.reuse_if_untouched(o)

    @process.register(ufl.classes.Terminal)
    def _(self, t: ufl.classes.Terminal) -> ufl.classes.Terminal:
        """Apply to a terminal."""
        return t

    @process.register(ufl.classes.FormArgument)
    def _(self, o: ufl.classes.FormArgument) -> ufl.classes.Expr:
        """Apply to a form_argument."""
        # Represent 0-derivatives of form arguments on reference
        # element
        r = ReferenceValue(o)
        space = o.ufl_function_space()
        element = o.ufl_element()

        if r.ufl_shape != element.reference_value_shape:
            raise ValueError(
                "Expecting reference space expression with shape "
                f"'{element.reference_value_shape}', got '{r.ufl_shape}'"
            )
        f = element.pullback.apply(r)
        if f.ufl_shape != space.value_shape:
            raise ValueError(
                f"Expecting pulled back expression with shape '{space.value_shape}', "
                f"got '{f.ufl_shape}'"
            )

        assert f.ufl_shape == o.ufl_shape
        return f


def apply_function_pullbacks(expr):
    """Change representation of coefficients and arguments in an expression.

    Applies Piola mappings where applicable and represents all
    form arguments in reference value.

    Args:
        expr: An Expression
    """
    return FunctionPullbackApplier()(expr)
