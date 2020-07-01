reticulate::use_python("/opt/anaconda3/bin/python")
library(reticulate)

pyo = import("pyomo.environ", convert = FALSE)

# had signal handling issues when using solver. Turned
# off signal handling based on this discussion thread
# https://github.com/PyUtilib/pyutilib/issues/31
pyulib = import("pyutilib.subprocess.GlobalData")
pyulib$DEFINE_SIGNAL_HANDLERS_DEFAULT = FALSE

source("operators.R")

M = pyo$ConcreteModel()
M$x = pyo$Var()
M$y = pyo$Var()

M$obj = pyo$Objective(expr = (M$x - 1)^2 + 100*(M$y - M$x^2)^2)

M$pprint()

solver = pyo$SolverFactory("ipopt")
res = solver$solve(M, tee = TRUE)

M$x();M$y()

