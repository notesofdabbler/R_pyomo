library(reticulate)
library(dplyr)
library(ggplot2)

pyo = import("pyomo.environ", convert = FALSE)
pyde = import("pyomo.dae", convert = FALSE)

source("operators.R")

m = pyo$ConcreteModel()

m$t = pyde$ContinuousSet(bounds = tuple(0, 1))
m$x1 = pyo$Var(m$t, bounds = tuple(0, 1))
m$x2 = pyo$Var(m$t, bounds = tuple(0, 1))
m$u = pyo$Var(m$t, initialize = 0)

m$x1dot = pyde$DerivativeVar(m$x1, wrt = m$t)
m$x2dot = pyde$DerivativeVar(m$x2, wrt = m$t)

m$obj = pyo$Objective(expr = m$x2[1], sense = pyo$minimize)

x1dot_rule = function(M, i) {
  if (i == 0) {
    return(py_get_attr(pyo$Constraint, "Skip"))
  } else {
    expr = M$x1dot[i] == M$u[i]
  }
  return(expr)
}
m$x1dot_cons = pyo$Constraint(m$t, rule = x1dot_rule)

x2dot_rule = function(M, i) {
  if (i == 0) {
    return(py_get_attr(pyo$Constraint, "Skip"))
  } else {
    expr = M$x2dot[i] == M$x1[i]^2 + M$u[i]^2
  }
  return(expr)
}

m$x2dot_cons = pyo$Constraint(m$t, rule = x2dot_rule)
m$init_x1 = pyo$Constraint(expr = m$x1[0] == 1.0)
m$init_x2 = pyo$Constraint(expr = m$x2[0] == 0.0)

discretizer = pyo$TransformationFactory('dae.collocation')
discretizer$apply_to(m,nfe=20L,ncp=3L,scheme='LAGRANGE-RADAU')
discretizer$reduce_collocation_points(m,var=m$u,ncp=1L,contset=m$t)

solver = pyo$SolverFactory("ipopt")
results = solver$solve(r_to_py(m), tee = TRUE)


get_pyo_resval = function(m, var) {
  outstr = 
    "
resval = []
for t in r.m.t:
    resval.append(r.m.x1[t]())
"
  outstr = gsub("m", m, outstr)
  outstr = gsub("x1", var, outstr)
  
  py_run_string(outstr)
  
  return(py$resval)
}

get_pyo_tval = function(m, var) {
  outstr = 
    "
tval = []
for t in r.m.t:
    tval.append(t)
"
  outstr = gsub("m", m, outstr)

  py_run_string(outstr)
  
  return(as.numeric(py$tval))
}

py_run_string("
tval = []
for t in r.m.t:
    tval.append(t) 
print(tval)
              ")


for (i in m$t){
  print(i)
}

resdf = tibble(t = get_pyo_tval("m"), 
               x1 = get_pyo_resval("m", "x1"), 
               x2 = get_pyo_resval("m", "x2"), 
               u = get_pyo_resval("m", "u"))

ggplot(resdf) + geom_line(aes(x = t, y = x1, color = "x1")) +
  geom_line(aes(x = t, y = x2, color = "x2")) +
  geom_line(aes(x = t, y = u, color = "u")) + theme_bw()
