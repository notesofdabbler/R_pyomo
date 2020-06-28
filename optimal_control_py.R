library(reticulate)

pyo = import("pyomo.environ")
pyde = import("pyomo.dae")

source("operators.R")

py_run_string("

import pyomo.environ as pyo
import pyomo.dae as pyde

m = pyo.ConcreteModel()

m.t = pyde.ContinuousSet(bounds = (0, 1))

m.x1 = pyo.Var(m.t, bounds = (0, 1))
m.x2 = pyo.Var(m.t, bounds = (0, 1))
m.u = pyo.Var(m.t, initialize = 0)

m.x1dot = pyde.DerivativeVar(m.x1, wrt = m.t)
m.x2dot = pyde.DerivativeVar(m.x2, wrt = m.t)

def _x1dot(M, i):
    if i == 0:
        return pyo.Constraint.Skip
    else:
        expr = M.x1dot[i] == m.u[i]
    return expr
m.x1dot_cons = pyo.Constraint(m.t, rule = _x1dot)

def _x2dot(M, i):
    if i == 0:
        return pyo.Constraint.Skip
    else:
        expr = M.x2dot[i] == m.x1[i]**2 + m.u[i]**2
    return expr
m.x2dot_cons = pyo.Constraint(m.t, rule = _x2dot)

m.init_x1 = pyo.Constraint(expr = m.x1[0] == 1)
m.init_x2 = pyo.Constraint(expr = m.x2[0] == 0)

#discretizer = pyo.TransformationFactory('dae.collocation')
#discretizer.apply_to(m,nfe=20,ncp=3,scheme='LAGRANGE-RADAU')
              ")

m2 = py_to_r(py$m)

discretizer = pyo$TransformationFactory('dae.collocation')
discretizer$apply_to(py$m,nfe=20L,ncp=3L,scheme='LAGRANGE-RADAU')

py$m$pprint()

discretizer = pyo$TransformationFactory('dae.collocation')
discretizer$apply_to(m2,nfe=20L,ncp=3L,scheme='LAGRANGE-RADAU')

m2$pprint()

solver = pyo$SolverFactory("ipopt")
results = solver$solve(m, tee = TRUE)

t = c()
x1 = c()
x2 = c()
u = c()
for (i in 1:61) {
  t[i] = m$t[i]
  x1[i] = m$x1[t[i]]()
  x2[i] = m$x2[t[i]]()
  u[i] = m$u[t[i]]()
}


