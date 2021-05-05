# -*- coding: utf-8 -*-
"""This module defines the Matrix class."""

# Copyright (C) 2021 India Marsden
#
# This file is part of UFL (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from ufl.form import BaseForm, FormSum, Form
from ufl.coefficient import Coefficient

# --- The Action class represents the adjoint of a numerical object that needs to be computed at compile time ---


class Action(BaseForm):
    """UFL base form type: respresents the action of an object on another"""

    __slots__ = (
        "_left",
        "_right",
        "_repr",
        "_arguments")
    _globalcount = 0

    def __getnewargs__(self):
        return (self._left, self._right)

    def __new__(cls, *args, **kw):
        assert(len(args) == 2)
        left = args[0]
        right = args[1]

        if isinstance(left, FormSum):
            # Adjoint distributes over sums on the LHS
            return FormSum(*[(Action(component, right), 1) for component in left.components()])
        if isinstance(right, FormSum):
            # Adjoint distributes over sums on the RHS
            return FormSum(*[(Action(left, component), 1) for component in right.components()])

        return super(Action, cls).__new__(cls)

    def __init__(self, left, right):
        BaseForm.__init__(self)

        self._left = left
        self._right = right

        if isinstance(right, Form):
            if self._left.arguments()[-1].ufl_function_space().dual() != self._right.arguments()[0].ufl_function_space():
                raise TypeError
        elif isinstance(right, Coefficient):
            if self._left.arguments()[-1].ufl_function_space() != self._right.ufl_function_space():
                raise TypeError
        else:
            raise TypeError

        self._repr = "Action(%s, %s)" % (repr(self._left), repr(self._right))

    def ufl_function_spaces(self):
        "Get the tuple of function spaces of the underlying form"
        if isinstance(self._right, Form):
            return self._left.ufl_function_spaces()[:-1] + self._right.ufl_function_spaces()[1:]
        elif isinstance(self._right, Coefficient):
            return self._left.ufl_function_spaces()[:-1]

    def _analyze_form_arguments(self):
        "Define arguments of a adjoint of a form as the reverse of the form arguments"
        if isinstance(self._right, Form):
            self._arguments = self._left.arguments()[:-1] + self._right.arguments()[1:]
        elif isinstance(self._right, Coefficient):
            self._arguments = self._left.arguments()[:-1]
        else:
            raise TypeError

    def __str__(self):
        return "Action(%s, %s)" % (repr(self._left), repr(self._right))

    def __repr__(self):
        return self._repr