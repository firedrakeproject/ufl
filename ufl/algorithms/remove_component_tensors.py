"""Remove component tensors.

This module contains classes and functions to remove component tensors.
"""
# Copyright (C) 2025 Pablo Brubeck
#
# This file is part of UFL (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import functools
from collections import defaultdict

import ufl.classes
from ufl.algorithms.map_integrands import map_integrand_dags
from ufl.classes import ComponentTensor, Index, MultiIndex, Zero
from ufl.corealg.dag_traverser import DAGTraverser
from ufl.corealg.map_dag import map_expr_dag
from ufl.index_combination_utils import unique_sorted_indices


class IndexReplacer(DAGTraverser):
    """Replace Indices."""

    def __init__(self, fimap: dict):
        """Initialise.

        Args:
           fimap: map for index replacements.

        """
        super().__init__()
        self.fimap = fimap

    @functools.singledispatchmethod
    def process(self, o: ufl.classes.Expr) -> ufl.classes.Expr:
        """Process ``o``."""
        return super().process(o)

    @process.register(ufl.classes.Expr)
    def _(self, o: ufl.classes.Expr) -> ufl.classes.Expr:
        return self.reuse_if_untouched(o)

    @process.register(ufl.classes.Zero)
    def _(self, o: ufl.classes.Zero) -> ufl.classes.Zero:
        """Handle Zero."""
        indices = tuple(map(Index, o.ufl_free_indices))
        if not any(i in self.fimap for i in indices):
            # Reuse if untouched
            return o

        fi = []
        for i, d in zip(indices, o.ufl_index_dimensions):
            j = self.fimap.get(i, i)
            if isinstance(j, Index):
                fi.append((j.count(), d))

        fi = unique_sorted_indices(sorted(fi))
        free_indices, index_dimensions = zip(*fi)

        return Zero(
            shape=o.ufl_shape,
            free_indices=free_indices,
            index_dimensions=index_dimensions,
        )

    @process.register(ufl.classes.MultiIndex)
    def _(self, o: ufl.classes.MultiIndex) -> ufl.classes.MultiIndex:
        """Handle MultiIndex."""
        if not any(i in self.fimap for i in o):
            # Reuse if untouched
            return o

        indices = tuple(self.fimap.get(i, i) for i in o)
        return MultiIndex(indices)


class IndexRemover(DAGTraverser):
    """Remove Indexed."""

    def __init__(self):
        """Initialise."""
        super().__init__()
        self.rules = {}

    @functools.singledispatchmethod
    def process(self, o: ufl.classes.Expr) -> ufl.classes.Expr:
        """Process ``o``."""
        return super().process(o)

    @process.register(ufl.classes.Expr)
    def _(self, o: ufl.classes.Expr) -> ufl.classes.Expr:
        return self.reuse_if_untouched(o)

    @process.register(ufl.classes.Indexed)
    @DAGTraverser.postorder
    def _(self, o: ufl.classes.Indexed, o1, i1) -> ufl.classes.Expr:
        """Simplify Indexed."""
        if isinstance(o1, ComponentTensor):
            # Simplify Indexed ComponentTensor
            o2, i2 = o1.ufl_operands
            # Replace outer indices
            rkey = (i2, i1)
            rule = self.rules.get(rkey)
            if rule is None:
                # NOTE: Replace with `fimap = dict(zip(i2, i1, strict=True))` when
                # Python>=3.10
                assert len(i2) == len(i1)
                fimap = dict(zip(i2, i1))
                rule = IndexReplacer(fimap)
                self.rules[rkey] = rule

            key = (IndexReplacer, *rkey)
            return self(o2)

        elif o.ufl_operands[0] is o1:
            # Reuse if untouched
            return o
        else:
            return o._ufl_expr_reconstruct_(o1, i1)


def remove_component_tensors(o):
    """Remove component tensors."""
    rule = IndexRemover()
    return rule(o)
