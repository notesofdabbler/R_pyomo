#'---
#'title: "Kinetics Parameter Estimation in  Pyomo using R"
#'author: "notesofdabbler"
#'---
#'
#' This example estimates reaction rates for the reactions $A \rightarrow B \rightarrow C$ 
#' based on experimental data on concentrations of $A$, $B$ and $C$ over time. This example
#' is taken from the book 
#' [Chemical Reactor Analysis and Design - Rawlings and Ekerdt](https://sites.engineering.ucsb.edu/~jbraw/chemreacfun/)
#' The parameter estimation problem is to find the rate constants $k1$, $k2$ and initial
#' concentrations of A, B, and C that minimize the sum of squares error between experimental
#' and predicted data.
#' $$
#' \min \;\; \sum_{i = 1}^N((ca(t_i) - ca^{meas}(t_i))^2 + 
#'                           (cb(t_i) - cb^{meas}(t_i))^2 +
#'                           (cc(t_i) - cc^{meas}(t_i))^2) \\
#'   \frac{dca}{dt} = -k_1ca \\
#'   \frac{dcb}{dt} = k_1ca - k_2cb \\
#'   \frac{dcc}{dt} = k_2cb \\
#'   ca[0] = ca0 \\
#'   cb[0] = cb0 \\
#'   cc[0] = cc0 \\
#'   0 \leq t \leq t_f
#' $$
#' where experimental data on concentrations is measured at $N$ time points $\{t_1, t_2, \ldots, t_N\}$
#+
# set-up python version to use
reticulate::use_python("/opt/anaconda3/bin/python")
library(reticulate)

# import pyomo
pyo = import("pyomo.environ", convert = FALSE)
pydae = import("pyomo.dae", convert = FALSE)
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
library(tidyr)
library(ggplot2)

# read exprimental data (this data was generated using simulation)
data_df = read.csv("python_codes/ABC_data.csv", stringsAsFactors = FALSE)
head(data_df)

# convert data for each component into a dictionary
data_df_py = r_to_py(data_df)$set_index('t')$to_dict()
ca_meas = data_df_py['ca']
cb_meas = data_df_py['cb']
cc_meas = data_df_py['cc']

# check the dictionary for component A
py_run_string("print(r.ca_meas)")

# create model object
M = pyo$ConcreteModel()

# initial concentration estimates
M$ca0 = pyo$Param(initialize = 1.0)
M$cb0 = pyo$Param(initialize = 0.0)
M$cc0 = pyo$Param(initialize = 0.0)

# time variable
M$t = pydae$ContinuousSet(bounds = tuple(0.0, 5.0), initialize = data_df$t)

# parameters: measurement times, measured concentrations ca, cb, cc
M$t_meas = pyo$Set(within = M$t, initialize = data_df$t)
M$ca_meas = pyo$Param(M$t_meas, initialize = ca_meas)
M$cb_meas = pyo$Param(M$t_meas, initialize = cb_meas)
M$cc_meas = pyo$Param(M$t_meas, initialize = cc_meas)

# variables: rate constants for the 2 reactions A->B, B->C
M$k1 = pyo$Var(initialize= 0.5, bounds = tuple(1e-4, 10))
M$k2 = pyo$Var(initialize = 3, bounds = tuple(1e-4, 10))

# variables: concentration of A, B, C over time
M$ca = pyo$Var(M$t, initialize = M$ca0, bounds = tuple(0.0, M$ca0))
M$cb = pyo$Var(M$t, initialize = M$cb0, bounds = tuple(0.0, M$ca0))
M$cc = pyo$Var(M$t, initialize = M$cc0, bounds = tuple(0.0, M$ca0))

# derivative variables of concentration (wrt time) for A, B, and C
M$ca_dot = pydae$DerivativeVar(M$ca, wrt = M$t)
M$cb_dot = pydae$DerivativeVar(M$cb, wrt = M$t)
M$cc_dot = pydae$DerivativeVar(M$cc, wrt = M$t)

# constraint: ca_dot = -k1*ca
dcarate_rule = function(m, i) {
   expr = m$ca_dot[i] == 0-m$k1 * m$ca[i]
   return(expr)
}
M$dcarate_cons = pyo$Constraint(M$t, rule = dcarate_rule)

# constraint: cb_dot = k1*ca - k2*cb
dcbrate_rule = function(m, t) {
  expr = m$cb_dot[t] == m$k1 * m$ca[t] - m$k2 * m$cb[t]
  return(expr)
}
M$dcbrate_cons = pyo$Constraint(M$t, rule = dcbrate_rule)

# constraint: cc_dot = k2*cb
dccrate_rule = function(m, t) {
  expr = m$cc_dot[t] == m$k2 * m$cb[t]
  return(expr)
}
M$dccrate_cons = pyo$Constraint(M$t, rule = dccrate_rule)

# objective: minimize sum of square error
# Note: I had to use tuple to index time point since otherwise it was giving an error
ssq_rule = function(m) {
  ssq = 0
  for (i in bi$list(m$t_meas)) {
    ssq = ssq + (m$ca[tuple(i)] - m$ca_meas[tuple(i)])**2 + (m$cb[tuple(i)] - m$cb_meas[tuple(i)])**2 +
      (m$cc[tuple(i)] - m$cc_meas[tuple(i)])**2
  }
  return(ssq)
}
M$ssq_obj = pyo$Objective(rule = ssq_rule, sense = pyo$minimize)

# apply collocation
disc = pyo$TransformationFactory('dae.collocation')
disc$apply_to(M, nfe=20L, ncp=2L)

# solve using IPOPT
solver = pyo$SolverFactory('ipopt')
solver$solve(M, logfile = "tmp.log")
cat(paste(readLines("tmp.log"), collapse = "\n"))

# Get the estimated parameters
k1_est = M$k1()
k2_est = M$k2()

print(paste0("k1 = ", round(py_to_r(k1_est), 2), "; k2 = ", round(py_to_r(k2_est), 2)))

# Create a tibble of simulated concentration profiles at the 
# estimated rate constants k1_est and k2_est
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

# plot experimental vs simulated data
p = ggplot()
p = p + geom_line(data = c_sim_df_l, aes(x = t_sim, y = conc, color = comp))
p = p + geom_point(data = data_df_l, aes(x = t, y = conc, color = comp))
p = p + labs(x = "time", y = "concentration")
p = p + theme_bw()
p

#'### Session Info
sessionInfo()
# pyomo version used is 5.7