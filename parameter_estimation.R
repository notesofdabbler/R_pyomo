reticulate::use_python("/Users/shanki/anaconda3/bin/python")
library(reticulate)
library(dplyr)
library(tidyr)
library(ggplot2)

source("operators.R")

pyo = import("pyomo.environ", convert = FALSE)
pydae = import("pyomo.dae", convert = FALSE)
bi = import_builtins()

data_df = read.csv("python_codes/ABC_data.csv", stringsAsFactors = FALSE)

data_df_py = r_to_py(data_df)$set_index('t')$to_dict()
ca_meas = data_df_py['ca']
cb_meas = data_df_py['cb']
cc_meas = data_df_py['cc']

py_run_string("print(r.ca_meas)")

M = pyo$ConcreteModel()

M$ca0 = pyo$Param(initialize = 1.0)
M$cb0 = pyo$Param(initialize = 0.0)
M$cc0 = pyo$Param(initialize = 0.0)

M$t = pydae$ContinuousSet(bounds = tuple(0.0, 5.0), initialize = data_df$t)
M$t_meas = pyo$Set(within = M$t, initialize = data_df$t)
M$ca_meas = pyo$Param(M$t_meas, initialize = ca_meas)
M$cb_meas = pyo$Param(M$t_meas, initialize = cb_meas)
M$cc_meas = pyo$Param(M$t_meas, initialize = cc_meas)

M$k1 = pyo$Var(initialize= 0.5, bounds = tuple(1e-4, 10))
M$k2 = pyo$Var(initialize = 3, bounds = tuple(1e-4, 10))

M$ca = pyo$Var(M$t, initialize = M$ca0, bounds = tuple(0.0, M$ca0))
M$cb = pyo$Var(M$t, initialize = M$cb0, bounds = tuple(0.0, M$ca0))
M$cc = pyo$Var(M$t, initialize = M$cc0, bounds = tuple(0.0, M$ca0))

M$ca_dot = pydae$DerivativeVar(M$ca, wrt = M$t)
M$cb_dot = pydae$DerivativeVar(M$cb, wrt = M$t)
M$cc_dot = pydae$DerivativeVar(M$cc, wrt = M$t)


dcarate_rule = function(m, i) {
   expr = m$ca_dot[i] == 0-m$k1 * m$ca[i]
   return(expr)
}
M$dcarate_cons = pyo$Constraint(M$t, rule = dcarate_rule)

dcbrate_rule = function(m, t) {
  expr = m$cb_dot[t] == m$k1 * m$ca[t] - m$k2 * m$cb[t]
  return(expr)
}
M$dcbrate_cons = pyo$Constraint(M$t, rule = dcbrate_rule)

dccrate_rule = function(m, t) {
  expr = m$cc_dot[t] == m$k2 * m$cb[t]
  return(expr)
}
M$dccrate_cons = pyo$Constraint(M$t, rule = dccrate_rule)

ssq_rule = function(m) {
  ssq = 0
  for (i in bi$list(m$t_meas)) {
    ssq = ssq + (m$ca[tuple(i)] - m$ca_meas[tuple(i)])**2 + (m$cb[tuple(i)] - m$cb_meas[tuple(i)])**2 +
      (m$cc[tuple(i)] - m$cc_meas[tuple(i)])**2
  }
  return(ssq)
}
M$ssq_obj = pyo$Objective(rule = ssq_rule, sense = pyo$minimize)

disc = pyo$TransformationFactory('dae.collocation')
disc$apply_to(M, nfe=20L, ncp=2L)

solver = pyo$SolverFactory('ipopt')
solver$solve(M, tee = TRUE)

k1_est = M$k1()
k2_est = M$k2()

t_sim = c()
ca_sim = c()
cb_sim = c()
cc_sim = c()

tlist = as.numeric(bi$list(M$t))
for(i in seq_along(tlist)) {
  t_sim[i] = tlist[i]
  ca_sim[i] = py_to_r(M$ca[tuple(tlist[i])]())
  cb_sim[i] = py_to_r(M$cb[tuple(tlist[i])]())
  cc_sim[i] = py_to_r(M$cc[tuple(tlist[i])]())
}

c_sim_df = tibble(t_sim, ca_sim, cb_sim, cc_sim)
c_sim_df = c_sim_df %>% rename(ca = ca_sim, cb = cb_sim, cc = cc_sim)
c_sim_df_l = pivot_longer(c_sim_df, ca:cc, names_to = "comp", values_to = "conc")
data_df_l = pivot_longer(data_df, ca:cc, names_to = "comp", values_to = "conc") 

p = ggplot()
p = p + geom_line(data = c_sim_df_l, aes(x = t_sim, y = conc, color = comp))
p = p + geom_point(data = data_df_l, aes(x = t, y = conc, color = comp))
p = p + theme_bw()
p
