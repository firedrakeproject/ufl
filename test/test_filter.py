import pytest
from ufl import *
from ufl.algorithms.formsplitter import extract_blocks
from ufl.algorithms import compute_form_data
from ufl.differentiation import ReferenceGrad


def test_filter_one_form():
    cell = triangle
    domain = Mesh(cell)
    element = VectorElement("Lagrange", cell, 1)
    V = FunctionSpace(domain, element)

    v = TestFunction(V)
    c = Coefficient(V)

    t_domain = TopologicalMesh(cell)
    t_V = TopologicalFunctionSpace(t_domain, element)
    fltr = Filter(t_V)

    filtered_v = Filtered(v, fltr)
    form = inner(c, grad(filtered_v[1])) * dx
    # form = inner(c, grad(filtered_v)[1,:]) * dx

    fd = compute_form_data(
        form,
        do_apply_function_pullbacks=True,
        do_apply_integral_scaling=True,
        do_apply_geometry_lowering=True,
        preserve_geometry_types=(),
        do_apply_restrictions=True,
        do_apply_filters=True,
        do_estimate_degrees=True,
        complex_mode=False
    )

    print()
    #itgd = fd0.integral_data[0].integrals[0].integrand()
    #x = SpatialCoordinate(domain)
    #F = ReferenceGrad(x)
    #_F = as_tensor([[2., 3.], [5., 7.]])
    #intgd = replace(intgd, {F: _F})
    #print(itgd.ufl_free_indices)
    print(fd)
    print(repr(fd.integral_data[0].integrals[0].integrand()))
    assert fd.num_topological_coefficients == 1
    assert fd.reduced_topological_coefficients[0] is fltr
    assert fd.original_topological_coefficient_positions == [0, ]
    assert fd.integral_data[0].integral_topological_coefficients == set((fltr, ))
    assert fd.integral_data[0].enabled_topological_coefficients == [True, ]
    """
    mesh = Mesh(VectorElement(FiniteElement('Lagrange', triangle, 1), dim=2), 3)
    FS = FunctionSpace(mesh, FiniteElement('Lagrange', triangle, 1))
    F = ReferenceGrad(SpatialCoordinate(mesh))
    F00 = Indexed(F, MultiIndex((FixedIndex(0), FixedIndex(0))))
    F01 = F01
    F10 = F10
    F11 = F11
    detF = Abs(Sum(Product(F00, F11), Product(IntValue(-1), Product(F01, F10))))

    Product(
      Product(QuadratureWeight(mesh), detF), 
      IndexSum(
        Product(
          Indexed(
            ComponentTensor(
              IndexSum(
                Product(
                  Indexed(
                    ComponentTensor(
                      Division(
                        Indexed(ListTensor(ListTensor(F11, Product(IntValue(-1), F01)), ListTensor(Product(IntValue(-1), F10), F00)), MultiIndex((Index(104), Index(105)))), 
                        Sum(Product(F00, F11), Product(IntValue(-1), Product(F01, F10)))
                      ), MultiIndex((Index(104), Index(105)))
                    ), MultiIndex((Index(103), Index(102)))
                  ), 
                  Indexed(
                    ReferenceGrad(
                      Filtered(
                        ReferenceValue(Argument(FunctionSpace(mesh, VectorElement(FiniteElement('Lagrange', triangle, 1), dim=2)), 0, None)), 
                        Filter(TopologicalFunctionSpace(TopologicalMesh(triangle, 0), VectorElement(FiniteElement('Lagrange', triangle, 1), dim=2)), 0)
                      )
                    ), MultiIndex((Index(101), Index(103)))
                  )
                ), MultiIndex((Index(103),))
              ), MultiIndex((Index(101), Index(102)))
            ), MultiIndex((FixedIndex(1), Index(99)))
          ), 
          Indexed(ReferenceValue(Coefficient(FunctionSpace(mesh, VectorElement(FiniteElement('Lagrange', triangle, 1), dim=2)), 53)), MultiIndex((Index(99),)))
        ), MultiIndex((Index(99),))
      )
    )
    """
