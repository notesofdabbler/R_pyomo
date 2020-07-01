#'---
#'title: "Minimize Rosenbrock Function with Pyomo Using R"
#'author: "notesofdabbler"
#'---
#+
#  Minimization of Rosenbrock function
#  taken from:
#  https://github.com/Pyomo/pyomo/blob/master/examples/doc/pyomobook/nonlinear-ch/rosen/rosenbrock.py
#  Jupyter notebook based on above example:
#  https://github.com/notesofdabbler/R_pyomo/blob/master/python_codes/rosen.ipynb
#

# set-up python version to use
reticulate::use_python("/opt/anaconda3/bin/python")
library(reticulate)

# import pyomo
pyo = import("pyomo.environ", convert = FALSE)

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

# create model object and decision variables
M = pyo$ConcreteModel()
M$x = pyo$Var(initialize = 1.5)
M$y = pyo$Var(initialize = 1.5)

# model objective to minimize
M$obj = pyo$Objective(expr = (M$x - 1)^2 + 100*(M$y - M$x^2)^2)

# print model
cat(mdl_print("M"))

# solve using ipopt
solver = pyo$SolverFactory("ipopt")
res = solver$solve(M, logfile = "tmp.log")
cat(paste(readLines("tmp.log"), collapse = "\n"))

# show the solution
M$x();M$y()

#'### Session Info
sessionInfo()
# pyomo version used is 5.7