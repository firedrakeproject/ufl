# -*- coding: utf-8 -*-
"""This module contains the apply_filters algorithm which propagates filters in a form towards the terminals."""

# Copyright (C) 2008-2016 Martin Sandve Alnæs
#
# This file is part of UFL (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later


from ufl.log import error
from ufl.classes import Argument, Filtered
from ufl.corealg.multifunction import MultiFunction
from ufl.corealg.map_dag import map_expr_dag
from ufl.algorithms.map_integrands import map_integrand_dags
from ufl.measure import integral_type_to_measure_name


class FilterRuleset(MultiFunction):
    def __init__(self, fltr):
        MultiFunction.__init__(self)
        self._fltr = fltr

    def terminal(self, o):
        return o

    expr = MultiFunction.reuse_if_untouched

    def reference_value(self, o):
        "Must act directly on reference value of argument objects."
        f, = o.ufl_operands
        assert f._ufl_is_terminal_
        assert isinstance(f, Argument)
        return Filtered(o, self._fltr)


class FilterRuleDispatcher(MultiFunction):
    def __init__(self):
        MultiFunction.__init__(self)

    def terminal(self, o):
        return o

    expr = MultiFunction.reuse_if_untouched

    def filtered(self, o, A, fltr):
        rules = FilterRuleset(fltr)
        return map_expr_dag(rules, A)


def apply_filters(expression):
    "Propagate filter nodes to wrap reference value argument directly."
    rules = FilterRuleDispatcher()
    return map_integrand_dags(rules, expression)
