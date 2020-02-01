#!/usr/bin/env py.test
# -*- coding: utf-8 -*-

from ufl import *

import pytest

from ufl import *
from ufl.domain import default_domain
from ufl.algorithms import compute_form_data
from ufl.algorithms.formsplitter import extract_blocks 


def test_mixed_element():
    
    gdim = 3

    cell3 = Cell("tetrahedron", geometric_dimension=gdim)
    cell2 = Cell("triangle", geometric_dimension=gdim)

    element3 = FiniteElement("BDM", cell3, 1)
    element2 = FiniteElement("Lagrange", cell2, 1)

    # Define a mixed element on a mixed cell
    mixed_element =  MixedElement(element3, element2, mixed=True)

    # Define a mixed domain
    domain = (Mesh(cell3),
              Mesh(cell2))

    V = FunctionSpace(domain, mixed_element)

    v = TestFunction(V)
    u = TrialFunction(V)
    f = Coefficient(V)


    assert(mixed_element.cell() == (cell3, cell2))
    assert(V.split() == (V[0], V[1]))
    assert(u.split() == (u[0], u[1]))
    assert(v.split() == (v[0], v[1]))
    assert(f.split() == (f[0], f[1]))

    n = FacetNormal(cell3)
    
    ds0 = Measure("ds", domain[0])
    a10 = dot(u[0], n) * v[1] * ds0
    #dx1 = Measure("dx", domain[1])
    #a10 = dot(u[0], n) * v[1] * dx1

    fd = compute_form_data(a10,
                           do_apply_function_pullbacks=True,
                           #do_apply_integral_scaling=True,
                           do_apply_geometry_lowering=True,
                           #preserve_geometry_types=(CellVolume, FacetArea),
                           #do_apply_restrictions=True,
                           #do_estimate_degrees=True,
                           #complex_mode=False
                          )


    for integral_data in fd.integral_data:
        for integral in integral_data.integrals:
            print("\nprinting a integral...\n")
            print(repr(integral))

"""    
from ufl.classes import (ReferenceValue,
                         Jacobian, JacobianInverse, JacobianDeterminant,
                         Index)



u0, u1 = TrialFunctions(W)
v0, v1 = TestFunctions(W)


#form = inner(grad(u0), grad(v0))*dx



"""

"""
def test_mixed_functionspace(self):
    # Domains
    domain_3d = default_domain(tetrahedron)
    domain_2d = default_domain(triangle)
    domain_1d = default_domain(interval)
    # Finite elements
    f_1d = FiniteElement("CG", interval, 1)
    f_2d = FiniteElement("CG", triangle, 1)
    f_3d = FiniteElement("CG", tetrahedron, 1)
    # Function spaces
    V_3d = FunctionSpace(domain_3d, f_3d)
    V_2d = FunctionSpace(domain_2d, f_2d)
    V_1d = FunctionSpace(domain_1d, f_1d)

    # MixedFunctionSpace = V_3d x V_2d x V_1d
    V = MixedFunctionSpace(V_3d, V_2d, V_1d)
    # Check sub spaces
    assert( V.num_sub_spaces() == 3 )
    assert( V.ufl_sub_space(0) == V_3d )
    assert( V.ufl_sub_space(1) == V_2d )
    assert( V.ufl_sub_space(2) == V_1d )

    # Arguments from MixedFunctionSpace
    (u_3d, u_2d, u_1d) = TrialFunctions(V)
    (v_3d, v_2d, v_1d) = TestFunctions(V)
    
    # Measures
    dx3 = Measure("dx", domain=V_3d)
    dx2 = Measure("dx", domain=V_2d)
    dx1 = Measure("dx", domain=V_1d)
    
    # Mixed variational form
    # LHS
    a_11 = u_1d*v_1d*dx1
    a_22 = u_2d*v_2d*dx2
    a_33 = u_3d*v_3d*dx3
    a_21 = u_2d*v_1d*dx1
    a_12 = u_1d*v_2d*dx1
    a_32 = u_3d*v_2d*dx2
    a_23 = u_2d*v_3d*dx2
    a_31 = u_3d*v_1d*dx1
    a_13 = u_1d*v_3d*dx1
    a = a_11 + a_22 + a_33 + a_21 + a_12 + a_32 + a_23 + a_31 + a_13
    # RHS
    f_1 = v_1d*dx1
    f_2 = v_2d*dx2
    f_3 = v_3d*dx3
    f = f_1 + f_2 + f_3

    # Check extract_block algorithm
    # LHS
    assert ( extract_blocks(a,0,0) == a_33 )
    assert ( extract_blocks(a,0,1) == a_23 )
    assert ( extract_blocks(a,0,2) == a_13 )
    assert ( extract_blocks(a,1,0) == a_32 )
    assert ( extract_blocks(a,1,1) == a_22 )
    assert ( extract_blocks(a,1,2) == a_12 )
    assert ( extract_blocks(a,2,0) == a_31 )
    assert ( extract_blocks(a,2,1) == a_21 )
    assert ( extract_blocks(a,2,2) == a_11 )
    # RHS
    assert ( extract_blocks(f,0) == f_3 )
    assert ( extract_blocks(f,1) == f_2 )
    assert ( extract_blocks(f,2) == f_1 )
"""
