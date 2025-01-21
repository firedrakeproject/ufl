import pytest

from ufl import (
    Coefficient,
    FacetNormal,
    FunctionSpace,
    Mesh,
    SpatialCoordinate,
    TestFunction,
    TrialFunction,
    adjoint,
    as_tensor,
    cofac,
    conditional,
    conj,
    derivative,
    ds,
    dx,
    grad,
    inner,
    tetrahedron,
)
from ufl.algorithms.check_arities import ArityMismatch
from ufl.algorithms.compute_form_data import compute_form_data
from ufl.finiteelement import FiniteElement
from ufl.pullback import identity_pullback
from ufl.sobolevspace import H1


def test_check_arities():
    # Code from bitbucket issue #49
    cell = tetrahedron
    D = Mesh(FiniteElement("Lagrange", cell, 1, (3,), identity_pullback, H1))
    V = FunctionSpace(D, FiniteElement("Lagrange", cell, 2, (3,), identity_pullback, H1))
    dv = TestFunction(V)
    du = TrialFunction(V)

    X = SpatialCoordinate(D)
    N = FacetNormal(D)

    u = Coefficient(V)
    x = X + u
    F = grad(x)
    n = cofac(F) * N

    M = inner(x, n) * ds
    L = derivative(M, u, dv)
    a = derivative(L, u, du)

    compute_form_data(M)
    compute_form_data(L)
    compute_form_data(a)


def test_complex_arities():
    cell = tetrahedron
    D = Mesh(FiniteElement("Lagrange", cell, 1, (3,), identity_pullback, H1))
    V = FunctionSpace(D, FiniteElement("Lagrange", cell, 2, (3,), identity_pullback, H1))
    v = TestFunction(V)
    u = TrialFunction(V)

    # Valid form.
    F = inner(u, v) * dx
    compute_form_data(F, complex_mode=True)
    # Check that adjoint conjugates correctly
    compute_form_data(adjoint(F), complex_mode=True)

    with pytest.raises(ArityMismatch):
        compute_form_data(inner(v, u) * dx, complex_mode=True)

    with pytest.raises(ArityMismatch):
        compute_form_data(inner(conj(v), u) * dx, complex_mode=True)


def test_product_arity():
    cell = tetrahedron
    D = Mesh(FiniteElement("Lagrange", cell, 1, (3,), identity_pullback, H1))
    V = FunctionSpace(D, FiniteElement("Lagrange", cell, 2, (3,), identity_pullback, H1))
    v = TestFunction(V)
    u = TrialFunction(V)

    with pytest.raises(ArityMismatch):
        F = inner(u, u) * dx
        compute_form_data(F, complex_mode=True)

    with pytest.raises(ArityMismatch):
        L = inner(v, v) * dx
        compute_form_data(L, complex_mode=False)


def test_zero_simplify_arity():
    cell = tetrahedron
    D = Mesh(FiniteElement("Lagrange", cell, 1, (3,), identity_pullback, H1))
    V = FunctionSpace(D, FiniteElement("Lagrange", cell, 2, (), identity_pullback, H1))
    v = TestFunction(V)
    u = Coefficient(V)

    nonzero = 1
    with pytest.raises(ArityMismatch):
        F = inner(u, v + nonzero) * dx
        compute_form_data(F)

    zero = as_tensor([0, u])[0]
    F = inner(u, v + zero) * dx
    compute_form_data(F)

    zero = conditional(u < 0, 0, 0)
    F = inner(u, v + zero) * dx
    compute_form_data(F)

    zero = conditional(u < 0, 0, conditional(u == 0, 0, 0))
    F = inner(u, v + zero) * dx
    compute_form_data(F)
