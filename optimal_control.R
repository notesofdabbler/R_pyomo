#'---
#'title: "Optimal Control Example in  Pyomo using R"
#'author: "notesofdabbler"
#'---
#+
#  Optimal control example taken from:
#  https://github.com/Pyomo/pyomo/blob/master/examples/dae/Optimal_Control.py
#  https://github.com/Pyomo/pyomo/blob/master/examples/dae/run_Optimal_Control.py
#
#
#	min X2(tf)
#	s.t.	X1_dot = u			
#		    X2_dot = X1^2 + u^2		
#		    tf = 1
#       X1(0) = 1, X2(0) = 0
#
# set-up python version to use
reticulate::use_python("/opt/anaconda3/bin/python")
library(reticulate)

# import pyomo
pyo = import("pyomo.environ", convert = FALSE)
pyde = import("pyomo.dae", convert = FALSE)
bi = import_builtins()

# had signal handling issues when using solver. Turned
# off signal handling based on this discussion thread
# https://github.com/PyUtilib/pyutilib/issues/31
pyulib = import("pyutilib.subprocess.GlobalData")
pyulib$DEFINE_SIGNAL_HANDLERS_DEFAULT = FALSE

# printing model to a file, reading it back into R to print in R markdown
mdl_print = function(m) {
  printstr = "
with open('tmp_model.txt', 'w') as f:
  r.mdl.pprint(ostream=f)
f.close()
              "
  printstr = gsub("mdl", m, printstr)
  py_run_string(printstr)
  printstr_R = paste(readLines("tmp_model.txt"), collapse = "\n")
  return(printstr_R)
}

# get operators
source("operators.R")

# load other libraries
library(dplyr)
library(ggplot2)

# create model object
M = pyo$ConcreteModel()

# set time variable
M$t = pyde$ContinuousSet(bounds = tuple(0, 1))

# set state variables x1, x2
M$x1 = pyo$Var(M$t, bounds = tuple(0, 1))
M$x2 = pyo$Var(M$t, bounds = tuple(0, 1))
# set control variable u
M$u = pyo$Var(M$t, initialize = 0)

# set variables denoting derivatives of state variables x1dot, x2dot
M$x1dot = pyde$DerivativeVar(M$x1, wrt = M$t)
M$x2dot = pyde$DerivativeVar(M$x2, wrt = M$t)

# objective: minimize x2 at tf = 1
M$obj = pyo$Objective(expr = M$x2[1], sense = pyo$minimize)

# define constraint x1dot = u
x1dot_rule = function(m, i) {
  if (i == 0) {
    return(py_get_attr(pyo$Constraint, "Skip"))
  } else {
    expr = m$x1dot[i] == m$u[i]
  }
  return(expr)
}
M$x1dot_cons = pyo$Constraint(M$t, rule = x1dot_rule)

# define constraint x2dot = x1^2 + u^2
x2dot_rule = function(m, i) {
  if (i == 0) {
    return(py_get_attr(pyo$Constraint, "Skip"))
  } else {
    expr = m$x2dot[i] == m$x1[i]^2 + m$u[i]^2
  }
  return(expr)
}
M$x2dot_cons = pyo$Constraint(M$t, rule = x2dot_rule)

# set initial conditions: x1(0) = 1, x2(0) = 0
M$init_x1 = pyo$Constraint(expr = M$x1[0] == 1.0)
M$init_x2 = pyo$Constraint(expr = M$x2[0] == 0.0)

# apply collocation
discretizer = pyo$TransformationFactory('dae.collocation')
discretizer$apply_to(M,nfe=20L,ncp=3L,scheme='LAGRANGE-RADAU')
discretizer$reduce_collocation_points(M,var=M$u,ncp=1L,contset=M$t)

solver = pyo$SolverFactory("ipopt")
results = solver$solve(M, logfile = "tmp.log")
cat(paste(readLines("tmp.log"), collapse = "\n"))

# Get the results
#
# NOTE: 
# the time points are in the list M$t but it needs to accessed using bi$list(M$t)
# state and control variables are index by elements of M$t. 
# However accessing say x1 at time point i using M$x[i]() doesn't work.
# What seems to work is using M$x[tuple(i)]()
#
res_L = list()
k = 1 
for(i in bi$list(M$t)) {
  res_L[[k]] = list(t = i, x1 = py_to_r(M$x1[tuple(i)]()),
                    x2 = py_to_r(M$x2[tuple(i)]()),
                    u = py_to_r(M$u[tuple(i)]()))
  k = k + 1
}
res_df = bind_rows(res_L)

# plot the state and control variables
p = ggplot(res_df) + geom_line(aes(x = t, y = x1, color = "x1")) 
p = p + geom_line(aes(x = t, y = x2, color = "x2")) 
p = p + geom_line(aes(x = t, y = u, color = "u"))
p = p + labs(x = "time", y = "")
p = p + theme_bw()
p

#'### Session Info
sessionInfo()
# pyomo version used is 5.7
