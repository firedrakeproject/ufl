# -*- coding: utf-8 -*-
"""This module defines the Filter class."""

# Copyright (C) 2008-2016 Martin Sandve Alnæs
#
# This file is part of UFL (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from ufl.log import error
from ufl.core.ufl_type import ufl_type
from ufl.core.terminal import Terminal
from ufl.finiteelement import FiniteElementBase
from ufl.domain import default_domain
from ufl.functionspace import AbstractFunctionSpace, FunctionSpace, MixedFunctionSpace
from ufl.split_functions import split
from ufl.utils.counted import counted_init

# --- The Filter class represents a filter applied to arguments ---


@ufl_type()
class Filter(Terminal):
    """UFL terminal type: Representation of a filter."""
    __slots__ = ("_count", "_ufl_function_space", "_repr")
    _globalcount = 0

    def __init__(self, function_space, count=None):
        Terminal.__init__(self)
        counted_init(self, count, Filter)


        if isinstance(function_space, FiniteElementBase):
            element = function_space
            domain = default_domain(element.cell())
            function_space = FunctionSpace(domain, element)
        elif not isinstance(function_space, AbstractFunctionSpace):
            error("Expecting a FunctionSpace or FiniteElement.")

        self._ufl_function_space = function_space

        self._repr = "Filter(%s, %s)" % (
            repr(self._ufl_function_space), repr(self._count))

    def count(self):
        return self._count

    def ufl_function_space(self):
        "Get the function space of this coefficient."
        return self._ufl_function_space

    def ufl_domain(self):
        "Shortcut to get the domain of the function space of this filter."
        return self._ufl_function_space.ufl_domain()

    def ufl_element(self):
        "Shortcut to get the finite element of the function space of this filter."
        return self._ufl_function_space.ufl_element()

    @property
    def ufl_shape(self):
        "This shall not be used."
        error("Filter has no shape (it is not a tensor expression).")

    @property
    def ufl_free_indices(self):
        "This shall not be used."
        error("Filter has no free indices (it is not a tensor expression).")

    @property
    def ufl_index_dimensions(self):
        "This shall not be used."
        error("Filter has no free indices (it is not a tensor expression).")

    def is_cellwise_constant(self):
        "Always True."
        return True

    def ufl_domains(self):
        "Return tuple of domains related to this terminal object."
        return ()

    def _ufl_signature_data_(self, renumbering):
        "Signature data for filters depend on the global numbering of the form arguments, filters and domains."
        count = renumbering[self]
        fsdata = self._ufl_function_space._ufl_signature_data_(renumbering)
        return ("Filter", count, fsdata)

    def __str__(self):
        count = str(self._count)
        if len(count) == 1:
            return "filter_%s" % count
        else:
            return "filter_{%s}" % count

    def __repr__(self):
        return self._repr

    def __eq__(self, other):
        if not isinstance(other, Filter):
            return False
        if self is other:
            return True
        return (self._count == other._count and
                self._ufl_function_space == other._ufl_function_space)
