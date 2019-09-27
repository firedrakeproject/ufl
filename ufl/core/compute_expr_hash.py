# -*- coding: utf-8 -*-
"""Non-recursive traversal-based hash computation algorithm.

Fast iteration over nodes in an ``Expr`` DAG to compute
memoized hashes for all unique nodes.
"""

# Copyright (C) 2015 Martin Sandve Alnæs
#
# This file is part of UFL.
#
# UFL is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFL. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Massimiliano Leoni, 2016


def compute_expr_hash(expr):
    """Compute hashes of *expr* and all its nodes efficiently, without using Python recursion."""
    if expr._hash is not None:
        return expr._hash
    # Postorder traversal, can't use unique_post_traversal, since that
    # uses a set which requires that this hash is computed.
    lifo = [(expr, list(expr.ufl_operands))]
    while lifo:
        expr, deps = lifo[-1]
        for i, dep in enumerate(deps):
            if dep is not None and dep._hash is None:
                lifo.append((dep, list(dep.ufl_operands)))
                deps[i] = None
                break
        else:
            if expr._hash is None:
                expr._hash = expr._ufl_compute_hash_()
            lifo.pop()
    return expr._hash
