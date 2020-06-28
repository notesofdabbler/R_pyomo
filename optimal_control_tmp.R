library(reticulate)

pyo = import("pyomo.environ", convert = FALSE)
pyde = import("pyomo.dae", convert = FALSE)

source("operators.R")

py_run_string("
import pyomo.environ as pyo
m = pyo.ConcreteModel()
              ")
#m = pyo$ConcreteModel()

m = py$mpy

py$m$t = pyde$ContinuousSet(bounds = tuple(0, 1))
py$m$x1 = pyo$Var(py$m$t, bounds = tuple(0, 1))
py$m$x2 = pyo$Var(py$m$t, bounds = tuple(0, 1))
py$m$u = pyo$Var(py$m$t, initialize = 0)

py$m$x1dot = pyde$DerivativeVar(py$m$x1, wrt = py$m$t)
py$m$x2dot = pyde$DerivativeVar(py$m$x2, wrt = py$m$t)

py$m$obj = pyo$Objective(expr = py$m$x2[1], sense = pyo$minimize)

x1dot_rule = function(M, i) {
  if (i == 0) {
    return(py_get_attr(pyo$Constraint, "Skip"))
  } else {
    expr = M$x1dot[i] == M$u[i]
  }
  return(expr)
}
py$m$x1dot_cons = pyo$Constraint(py$m$t, rule = x1dot_rule)

x2dot_rule = function(M, i) {
  if (i == 0) {
    return(py_get_attr(pyo$Constraint, "Skip"))
  } else {
    expr = M$x2dot[i] == M$x1[i]^2 + M$u[i]^2
  }
  return(expr)
}
py$m$x2dot_cons = pyo$Constraint(py$m$t, rule = x2dot_rule)

py$m$init_x1 = pyo$Constraint(expr = py$m$x1[0] == 1.0)
py$m$init_x2 = pyo$Constraint(expr = py$m$x2[0] == 0.0)

m2 = r_to_py(m)

py_run_string("
print(sorted(r.m2.t))  
r.m2.x1dot_cons.pprint()
              ")

# Discretize model using Orthogonal Collocation
py_run_string("
import pyomo.environ as pyo
import copy

m3 = copy.deepcopy(r.m2)
m3.pprint()

# Discretize model using Backward Finite Difference method
discretizer = pyo.TransformationFactory('dae.finite_difference')
discretizer.apply_to(m3,nfe=20,scheme='BACKWARD')


#discretizer = pyo.TransformationFactory('dae.collocation')
#discretizer.apply_to(m3,nfe=20,ncp=3,scheme='LAGRANGE-RADAU')
#discretizer.reduce_collocation_points(m3,var=m3.u,ncp=1,contset=m3.t)

m3.pprint()
              ")



discretizer = pyo$TransformationFactory('dae.collocation')
discretizer$apply_to(py$m,nfe=20L,ncp=3L,scheme='LAGRANGE-RADAU')
discretizer$reduce_collocation_points(m2,var=m2$u,ncp=1L,contset=m2$t)



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


