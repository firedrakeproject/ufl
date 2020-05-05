# -*- coding: utf-8 -*-
"""This module defines the Transformed class."""

# Copyright (C) 2008-2016 Martin Sandve Alnæs
#
# This file is part of UFL (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from ufl.constantvalue import Zero
from ufl.core.ufl_type import ufl_type
from ufl.core.operator import Operator
from ufl.precedence import parstr


@ufl_type(num_ops=2, is_terminal_modifier=True, inherit_shape_from_operand=0, inherit_indices_from_operand=0)
class Transformed(Operator):
    __slots__ = (
        "ufl_shape",
        "ufl_free_indices",
        "ufl_index_dimensions",
    )

    def __new__(cls, expression, transform_op):
        if isinstance(expression, Zero):
            # Zero-simplify indexed Zero objects
            shape = expression.ufl_shape
            fi = expression.ufl_free_indices
            fid = expression.ufl_index_dimensions
            return Zero(shape=shape, free_indices=fi, index_dimensions=fid)
        else:
            return Operator.__new__(cls)

    def __init__(self, expression, transform_op):
        # Store operands
        Operator.__init__(self, (expression, transform_op))

        # Error checking
        #if not isinstance(expression, Argument):
        #    error("Expecting Argument instance, not %s." % ufl_err_str(expression))
        #if not isinstance(transform_op, TopologicalCoefficient):
        #    error("Expecting a TopologicalCoefficient instance, not %s." % ufl_err_str(transform_op))

        # Error checking
        #if len(shape) != len(multiindex):
        #    error("Invalid number of indices (%d) for tensor "
        #          "expression of rank %d:\n\t%s\n"
        #          % (len(multiindex), len(expression.ufl_shape), ufl_err_str(expression)))
        #if any(int(di) >= int(si)
        #       for si, di in zip(shape, multiindex)
        #       if isinstance(di, FixedIndex)):
        #    error("Fixed index out of range!")

        # Cache free index and dimensions
        #self.ufl_shape = expression.ufl_shape
        #self.ufl_free_indices = expression.ufl_free_indices
        #self.ufl_index_dimensions = expression.ufl_index_dimensions

    def ufl_element(self):
        "Shortcut to get the finite element of the function space of the operand."
        return self.ufl_operands[0].ufl_element()

    def __str__(self):
        return "%s[%s]" % (parstr(self.ufl_operands[0], self),
                           self.ufl_operands[1])
