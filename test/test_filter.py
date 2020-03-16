import pytest
from ufl import *
from ufl.algorithms.formsplitter import extract_blocks
from ufl.algorithms import compute_form_data as ufl_compute_form_data
from ufl.differentiation import ReferenceGrad


def test_filter_one_form():
    cell = triangle
    domain = Mesh(cell)
    element = VectorElement("Lagrange", cell, 1)
    V = FunctionSpace(domain, element)

    v = TestFunction(V)
    c = Coefficient(V)

    fltr = Filter(V)

    filtered_v = Filtered(v, fltr)
    form = inner(c, grad(filtered_v[1])) * dx
    # form = inner(c, grad(filtered_v)[1,:]) * dx

    fd = ufl_compute_form_data(
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
    assert fd.num_filters == 1
    assert fd.reduced_filters[0] is fltr
    assert fd.original_filter_positions == [0, ]
    assert fd.integral_data[0].integral_filters == set((fltr, ))
    assert fd.integral_data[0].enabled_filters == [True, ]
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
                        Indexed(
                          ListTensor(
                            ListTensor(F11, Product(IntValue(-1), F01)), 
                            ListTensor(Product(IntValue(-1), F10), F00)
                          ), MultiIndex((Index(104), Index(105)))
                        ), 
                        Sum(Product(F00, F11), Product(IntValue(-1), Product(F01, F10)))
                      ), MultiIndex((Index(104), Index(105)))
                    ), MultiIndex((Index(103), Index(102)))
                  ), 
                  Indexed(
                    ReferenceGrad(
                      Filtered(
                        ReferenceValue(Argument(FunctionSpace(mesh, VectorElement(FiniteElement('Lagrange', triangle, 1), dim=2)), 0, None)), 
                        Filter(FunctionSpace(mesh, VectorElement(FiniteElement('Lagrange', triangle, 1), dim=2)), 0)
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
