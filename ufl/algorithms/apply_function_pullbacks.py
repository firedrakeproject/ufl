"""Algorithm for replacing gradients in an expression."""

# Copyright (C) 2008-2016 Martin Sandve Alnæs
#
# This file is part of UFL (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from __future__ import annotations

from functools import singledispatchmethod

from ufl.algorithms.map_integrands import map_integrand_dags, map_integrands
from ufl.classes import Expr, Interpolate, ReferenceValue
from ufl.corealg.dag_traverser import DAGTraverser
from ufl.corealg.multifunction import MultiFunction, memoized_handler
from ufl.domain import AbstractDomain, extract_unique_domain
from ufl.finiteelement import AbstractFiniteElement
from ufl.form import BaseForm


class FunctionPullbackApplier(MultiFunction):
    """A pull back applier."""

    def __init__(self):
        """Initalise."""
        MultiFunction.__init__(self)

    expr = MultiFunction.reuse_if_untouched

    def terminal(self, t):
        """Apply to a terminal."""
        return t

    @memoized_handler
    def form_argument(self, o):
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


class InversePullbackApplier(DAGTraverser):
    """An inverse pull back applier.

    Args:
        element: The element whose pull back is inverted.
        domain: The domain to use if an expression carries none.
        compress: If True, ``result_cache`` will be used.
        visited_cache: cache of intermediate results; expr -> r = self.process(expr, ...).
        result_cache: cache of result objects for memory reuse, r -> r.

    """

    def __init__(
        self,
        element: AbstractFiniteElement,
        domain: AbstractDomain | None = None,
        compress: bool | None = True,
        visited_cache: dict[tuple, Expr | BaseForm] | None = None,
        result_cache: dict[Expr | BaseForm, Expr | BaseForm] | None = None,
    ) -> None:
        """Initialise."""
        super().__init__(compress=compress, visited_cache=visited_cache, result_cache=result_cache)
        self._element = element
        self._domain = domain

    @singledispatchmethod
    def process(self, o: Expr) -> Expr:
        """Map an expression onto the reference cell of ``self._element``.

        Args:
            o: An expression on a physical cell, whose shape must be the
                physical value shape of ``self._element``.

        Returns:
            The expression on the reference cell, with shape
            ``self._element.reference_value_shape``.

        """
        return super().process(o)

    @process.register(Expr)
    def _(self, o: Expr) -> Expr:
        """Handle Expr."""
        element = self._element
        mesh = extract_unique_domain(o) or self._domain
        if self._domain is not None and mesh != self._domain:
            raise NotImplementedError("Multiple domains not supported")
        physical_value_shape = element.pullback.physical_value_shape(element, mesh)
        if o.ufl_shape != physical_value_shape:
            raise ValueError(
                f"Expecting physical expression with shape '{physical_value_shape}', "
                f"got '{o.ufl_shape}'"
            )
        r = element.pullback.apply_inverse(o, mesh)
        if r.ufl_shape != element.reference_value_shape:
            raise ValueError(
                f"Expecting reference expression with shape "
                f"'{element.reference_value_shape}', got '{r.ufl_shape}'"
            )
        return r


class InterpolatePullbackApplier(DAGTraverser):
    """A pull back applier for interpolation."""

    @singledispatchmethod
    def process(self, o: Expr | BaseForm) -> Expr | BaseForm:
        """Process ``o``.

        Args:
            o: `Expr` or `BaseForm` to be processed.

        Returns:
            Processed `Expr` or `BaseForm`.

        """
        return super().process(o)

    @process.register(Expr)
    @process.register(BaseForm)
    def _(self, o: Expr | BaseForm) -> Expr | BaseForm:
        """Handle Expr and BaseForm."""
        return self.reuse_if_untouched(o)

    @process.register(Interpolate)
    @DAGTraverser.postorder
    def _(self, o: Interpolate, operand: Expr) -> Expr:
        """Evaluate an interpolation on the reference cell of its target element."""
        dual_arg, _ = o.argument_slots()
        element = o.ufl_element()
        domain = extract_unique_domain(operand) or dual_arg.ufl_function_space().ufl_domain()
        # Build the node here rather than reconstructing o: the mapped operand
        # no longer has the physical value shape that a subclass may check.
        r = Interpolate(apply_inverse_pullback(operand, element, domain), dual_arg)
        return element.pullback.apply(ReferenceValue(r), domain)


def apply_inverse_pullback(expr, element, domain=None):
    """Map a physical expression onto the reference cell of an element.

    Args:
        expr: An expression on a physical cell, whose shape must be the
            physical value shape of the element
        element: The element whose pull back is inverted
        domain: The domain to use if the expression carries none

    Returns:
        The expression on the reference cell, with shape
        ``element.reference_value_shape``
    """
    return InversePullbackApplier(element, domain=domain)(expr)


def apply_interpolate_pullbacks(expr):
    """Change the representation of the interpolations in an expression.

    An interpolation is evaluated on the reference cell of its target element,
    so its operand is mapped there and the result is pulled back to the
    physical cell for the expression that holds it.

    Args:
        expr: An Expr or Form

    Returns:
        The expression with its interpolations on their reference cells
    """
    return map_integrands(InterpolatePullbackApplier(), expr)


def apply_function_pullbacks(expr):
    """Change representation of coefficients and arguments in an expression.

    Applies Piola mappings where applicable and represents all
    form arguments in reference value.

    Args:
        expr: An Expression
    """
    return map_integrand_dags(FunctionPullbackApplier(), expr)
