# -*- coding: utf-8 -*-
"""Algorithm for spliting mixed coefficients."""

# Copyright (C) 2008-2016 Martin Sandve Alnæs and Anders Logg
#
# This file is part of UFL (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from ufl.log import error
from ufl.classes import CoefficientDerivative
from ufl.constantvalue import as_ufl
from ufl.core.terminal import FormArgument
from ufl.core.multiindex import FixedIndex
from ufl.corealg.multifunction import MultiFunction
from ufl.algorithms.map_integrands import map_integrand_dags
from ufl.algorithms.analysis import has_exact_type


class MixedCoefficientSplitRuleset(MultiFunction):
    def __init__(self):
        MultiFunction.__init__(self)

    def indexed(self, o, f, ii):
        if isinstance(f, FormArgument) and f.mixed():
            indices = ii.indices()
            assert len(indices) == 1, ""
            i = indices[0]
            assert isinstance(i, FixedIndex), ""
            i = int(i)
            return f.split()[i]
        else:
            return o

    expr = MultiFunction.reuse_if_untouched


def split_mixed_coefficients(expression):
    """Split mixed coefficients into component coefficients.

    @param expression:
        An Expr or Form.
    """
    rules = MixedCoefficientSplitRuleset()
    return map_integrand_dags(rules, expression)
