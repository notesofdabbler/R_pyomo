library(reticulate)

pyo = import("pyomo.environ", convert = FALSE)

source("operators.R")

M = pyo$ConcreteModel()
M$x = pyo$Var()
M$y = pyo$Var()

M$obj = pyo$Objective(expr = (M$x - 1)^2 + (M$y - M$x^2)^2)

M$pprint()

solver = pyo$SolverFactory("ipopt")
res = solver$solve(M, tee = TRUE)

M$x();M$y()

