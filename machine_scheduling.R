#
#  Machine scheduling problem:
#  https://nbviewer.jupyter.org/github/jckantor/ND-Pyomo-Cookbook/blob/master/notebooks/04.01-Machine-Bottleneck.ipynb
#
#

library(reticulate)
library(dplyr)
library(ggplot2)

pyo = import("pyomo.environ", convert = FALSE)
pygdp = import("pyomo.gdp", convert = FALSE)
bi = import_builtins()

source("operators.R")

jobs_df = tibble(job = c("A", "B", "C", "D", "E", "F", "G"),
                 release = c(2, 5, 4, 0, 0, 8, 9),
                 duration = c(5, 6, 8, 4, 2, 3, 2),
                 due = c(10, 21, 15, 10, 5, 15, 22))

release_L = list()
duration_L = list()
due_L = list()
for (i in 1: nrow(jobs_df)) {
  release_L[[jobs_df$job[i]]] = jobs_df$release[i]
  duration_L[[jobs_df$job[i]]] = jobs_df$duration[i]
  due_L[[jobs_df$job[i]]] = jobs_df$due[i]
}

tmax = sum(as.numeric(release_L)) + sum(as.numeric(duration_L))

M = pyo$ConcreteModel()

M$jobs = pyo$Set(initialize = jobs_df$job)

filt_rule = function(m, j, k) {
    return(j < k)
}
M$pairs = pyo$Set(initialize = M$jobs * M$jobs, filter = filt_rule)

M$release = pyo$Param(M$jobs, initialize = release_L)
M$duration = pyo$Param(M$jobs, initialize = duration_L)
M$due = pyo$Param(M$jobs, initialize = due_L)

M$start = pyo$Var(M$jobs, bounds = tuple(0, tmax))
M$pastdue = pyo$Var(M$jobs, bounds = tuple(0, tmax))
M$early = pyo$Var(M$jobs, bounds = tuple(0, tmax))

# ovlp_list_rule = function(b, j, k, flag) {
#   m2 = b$model()
#   if (flag == 0) {
#     b$c = pyo$Constraint(expr = m2$start[j] + m2$duration[j] <= m2$start[k])
#   } else {
#     b$c = pyo$Constraint(expr = m2$start[k] + m2$duration[k] <= m2$start[j])
#   }
# }
# M$ovlp_list_rule = pygdp$Disjunct(M$jobs, M$jobs, c(0, 1), rule = ovlp_list_rule)

ovlp_rule = function(m, j, k) {
  expr1 = m$start[j] + m$duration[j] <= m$start[k]
  expr2 = m$start[k] + m$duration[k] <= m$start[j]
  return(c(expr1, expr2))
}
M$ovlp_cons = pygdp$Disjunction(M$pairs, rule = ovlp_rule)


obj_rule = function(m) {
  expr = 0
  for (j in bi$list(m$jobs)) {
    expr = expr + m$pastdue[j]
  }
  return(expr)
}
M$obj = pyo$Objective(rule = obj_rule, sense = pyo$minimize)

pastdue_calc_rule = function(m, j) {
  expr = m$start[j] + m$duration[j] + m$early[j] == m$due[j] + m$pastdue[j]
  return(expr)
}
M$pastdue_calc_cons = pyo$Constraint(M$jobs, rule = pastdue_calc_rule)

release_rule = function(m, j) {
  expr = m$start[j] >= m$release[j]
  return(expr)
}
M$release_cons = pyo$Constraint(M$jobs, rule = release_rule)

pyo$TransformationFactory("gdp.chull")$apply_to(M)
res = pyo$SolverFactory('glpk')$solve(M)

resL = list()
for (i in bi$list(M$jobs)) {
  resL[[i]] = c(job = i, start = py_to_r(M$start[i]()), pastdue = py_to_r(M$pastdue[i]()))
}
res_df = bind_rows(resL)
