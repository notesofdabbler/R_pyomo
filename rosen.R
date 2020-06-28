library(reticulate)

"+.python.builtin.object" = function(a, b) {
  return(op$add(a, b))
}

"-.python.builtin.object" = function(a, b) {
  return(op$sub(a, b))
}

"*.python.builtin.object" = function(a, b) {
  return(op$mul(a, b))
}

"/.python.builtin.object" = function(a, b) {
  return(op$truediv(a, b))
}

"^.python.builtin.object" = function(a, b) {
  return(op$pow(a, b))
}

pyo = import("pyomo.environ")
op = import("operator")

M = pyo$ConcreteModel()
M$x = pyo$Var()
M$y = pyo$Var()

tmp = op$sub(M$x, 1)
tmp2 = op$pow(tmp, 2)
tmp3 = op$pow(op$sub(M$y, op$pow(M$x, 2)), 2)
tmp4 = op$add(tmp2, tmp3)

M$obj = pyo$Objective(expr = (M$x - 1)^2 + (M$y - M$x^2)^2)

M$pprint()

solver = pyo$SolverFactory("ipopt")
res = solver$solve(M, tee = TRUE)

M$x();M$y()

M1 = pyo$ConcreteModel()
M1$x = pyo$Var()
M1$y = pyo$Var()

cons = function(m) {
  tmp = op$sub(m$x, 1)
  tmp2 = op$pow(tmp, 2)
  tmp3 = op$pow(op$sub(m$y, op$pow(m$x, 2)), 2)
  tmp4 = op$add(tmp2, tmp3)
  return(tmp4)
}
M1$obj = pyo$Objective(rule = cons)

M1$pprint()

solver = pyo$SolverFactory("ipopt")
res = solver$solve(M1, tee = TRUE)
