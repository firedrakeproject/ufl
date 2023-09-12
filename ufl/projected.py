# -*- coding: utf-8 -*-
"""This module defines the Projected class."""

from ufl.constantvalue import Zero
from ufl.core.ufl_type import ufl_type
from ufl.core.operator import Operator
from ufl.subspace import Subspace
from ufl.precedence import parstr


# --- Indexed expression ---

@ufl_type(num_ops=2, is_terminal_modifier=True, inherit_shape_from_operand=0, inherit_indices_from_operand=0)
class Projected(Operator):
    __slots__ = (
        "ufl_shape",
        "ufl_free_indices",
        "ufl_index_dimensions",
    )

    def __new__(cls, expression, subspace):
        if isinstance(expression, Zero):
            # Zero-simplify indexed Zero objects
            shape = expression.ufl_shape
            fi = expression.ufl_free_indices
            fid = expression.ufl_index_dimensions
            return Zero(shape=shape, free_indices=fi, index_dimensions=fid)
        else:
            return Operator.__new__(cls)

    def __init__(self, expression, subspace):
        # Store operands
        Operator.__init__(self, (expression, subspace))
        # Error checking
        if not isinstance(expression, Expr):
            raise ValueError(f"Expecting Expr instance, not {ufl_err_str(expression)}.")
        if not isinstance(subspace, Subspace):
            raise ValueError(f"Expecting Subspace instance, not {ufl_err_str(subspace)}.")

    # def ufl_element(self):
    #    "Shortcut to get the finite element of the function space of the operand."
    #    return self.ufl_operands[0].ufl_element()

    def __str__(self):
        return "%s[%s]" % (parstr(self.ufl_operands[0], self),
                           self.ufl_operands[1])
