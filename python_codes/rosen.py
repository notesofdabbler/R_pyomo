#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example from
https://github.com/Pyomo/pyomo/blob/master/examples/pyomo/concrete/rosen.py
"""


#
# pyomo/core/expr/numvalue.py
# has functions for pyomo intrinsic operations
#
#

# rosen.py
from pyomo.environ import *
from pyomo.core.expr.current import *

M = ConcreteModel()
M.x = Var()
M.y = Var()
expr1 = (M.x-1)**2 + 100*(M.y-M.x**2)**2
tmp = M.x.__sub__(1)
tmp2 = pow(tmp, 2)
tmp3 = pow(M.y.__sub__(pow(M.x, 2)), 2)
tmp4 = tmp2.__add__(tmp3)


M.o  = Objective(
          expr=tmp4)

solver = SolverFactory('ipopt')
res = solver.solve(M, tee = True)

M.x()