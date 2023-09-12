# -*- coding: utf-8 -*-
"""This module defines the Subspace class."""

from ufl.core.ufl_type import ufl_type
from ufl.core.terminal import Terminal
from ufl.finiteelement import FiniteElementBase
from ufl.domain import default_domain
from ufl.functionspace import AbstractFunctionSpace, FunctionSpace, MixedFunctionSpace
from ufl.form import BaseForm
from ufl.split_functions import split
from ufl.utils.counted import Counted
from ufl.duals import is_primal, is_dual

# --- The Subspace class represents a subspace in a form ---


class BaseSubspace(Counted):
    """UFL terminal type: Parent Representatio/n of a subspace."""

    # Slots are disabled here because they cause trouble in PyDOLFIN
    # multiple inheritance pattern:
    # __slots__ = ("_count", "_ufl_function_space", "_repr", "_ufl_shape")
    _ufl_noslots_ = True
    __slots__ = ()
    _ufl_is_abstract_ = True

    def __getnewargs__(self):
        return (self._ufl_function_space, self._count)

    def __init__(self, function_space, count=None):
        Counted.__init__(self, count, Subspace)

        if isinstance(function_space, FiniteElementBase):
            # For legacy support for .ufl files using cells, we map
            # the cell to The Default Mesh
            element = function_space
            domain = default_domain(element.cell())
            function_space = FunctionSpace(domain, element)
        elif not isinstance(function_space, AbstractFunctionSpace):
            raise ValueError("Expecting a FunctionSpace or FiniteElement.")

        self._ufl_function_space = function_space
        ##self._ufl_shape = function_space.ufl_element().value_shape()

        self._repr = "BaseSubspace(%s, %s)" % (
            repr(self._ufl_function_space), repr(self._count))

    @property
    def ufl_shape(self):
        "This shall not be used."
        raise ValueError("Subspace has no shape (it is not a tensor expression).")

    @property
    def ufl_free_indices(self):
        "This shall not be used."
        raise ValueError("Subspace has no free indices (it is not a tensor expression).")

    @property
    def ufl_index_dimensions(self):
        "This shall not be used."
        raise ValueError("Subspace has no free indices (it is not a tensor expression).")

    #@property
    #def ufl_shape(self):
    #    "Return the associated UFL shape."
    #    return self._ufl_shape

    def ufl_function_space(self):
        "Get the function space of this subspace."
        return self._ufl_function_space

    def ufl_domain(self):
        "Shortcut to get the domain of the function space of this subspace."
        return self._ufl_function_space.ufl_domain()

    def ufl_element(self):
        "Shortcut to get the finite element of the function space of this subspace."
        return self._ufl_function_space.ufl_element()

    def is_cellwise_constant(self):
        "Always True. probably does not matter"
        return True

    #def is_cellwise_constant(self):
    #    "Return whether this expression is spatially constant over each cell."
    #    return self.ufl_element().is_cellwise_constant()

    def ufl_domains(self):
        "Return tuple of domains related to this terminal object."
        return self._ufl_function_space.ufl_domains()

    def _ufl_signature_data_(self, renumbering):
        "Signature data for form arguments depend on the global numbering of the form arguments and domains."
        count = renumbering[self]
        fsdata = self._ufl_function_space._ufl_signature_data_(renumbering)
        return ("Subspace", count, fsdata)

    def __str__(self):
        return f"s_{self._count}"

    def __repr__(self):
        return self._repr

    def __eq__(self, other):
        if not isinstance(other, BaseSubspace):
            return False
        if self is other:
            return True
        return self._count == other._count and self._ufl_function_space == other._ufl_function_space


@ufl_type()
class Subspace(Terminal, BaseSubspace]):
    """UFL terminal type: Representation of a subspace."""

    _ufl_noslots_ = True
    #_primal = True
    #_dual = False

    __getnewargs__ = BaseSubspace.__getnewargs__
    __str__ = BaseSubspace.__str__
    _ufl_signature_data_ = BaseSubspace._ufl_signature_data_

    def __init__(self, function_space, count=None):
        Terminal.__init__(self)
        BaseSubspace.__init__(self, function_space, count)

        self._repr = "Subspace(%s, %s)" % (
            repr(self._ufl_function_space), repr(self._count))

    def ufl_domains(self):
        return BaseSubspace.ufl_domains(self)

    def __eq__(self, other):
        if not isinstance(other, Subspace):
            return False
        if self is other:
            return True
        return self._count == other._count and self._ufl_function_space == other._ufl_function_space

    def __repr__(self):
        return self._repr
