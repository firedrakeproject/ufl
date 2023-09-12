#!/usr/bin/env py.test
# -*- coding: utf-8 -*-

"""
Test use of Projected.
"""

import pytest
from ufl import *
from ufl.algorithms import compute_form_data


def test_projected_one_form():
    cell = triangle
    element = VectorElement("Lagrange", cell, 1)
    domain = Mesh(cell)
    V = FunctionSpace(domain, element)
    V0 = Subspace(V)
    c = Coefficient(V)
    v = TestFunction(V)
    v0 = Projected(v, V0)
    form = inner(c, grad(v0[1])) * dx
    fd = compute_form_data(
        form,
        do_apply_function_pullbacks=True,
        do_apply_integral_scaling=True,
        do_apply_geometry_lowering=True,
        preserve_geometry_types=(),
        do_apply_restrictions=True,
        do_apply_projections=True,
        do_estimate_degrees=True,
        complex_mode=True
    )
    assert fd.num_subspaces == 1
    assert fd.reduced_subspaces[0] is V0
    assert fd.original_subspace_positions == [0, ]
    assert fd.integral_data[0].integral_subspaces == set((V0, ))
    assert fd.integral_data[0].enabled_subspaces == [True, ]
