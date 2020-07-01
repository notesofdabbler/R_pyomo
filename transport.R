#'---
#'title: "Solve Min Cost Transportation Problem with Pyomo using R"
#'author: "notesofdabbler"
#'---
#+
#  Min cost transportation problem
#  taken from:
#  https://nbviewer.jupyter.org/github/Pyomo/PyomoGallery/blob/master/transport/transport.ipynb
#
# set-up python version to use
reticulate::use_python("/opt/anaconda3/bin/python")
library(reticulate)

# import pyomo
pyo = import("pyomo.environ", convert = FALSE)
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

# shipping cost between plants and markets
costdf = tribble(
  ~plants, ~mkts, ~dtab,
  "seattle", "new-york", 2.5,
  "seattle", "chicago", 1.7,
  "seattle", "topeka", 1.8,
  "san-diego", "new-york", 2.5,
  "san-diego", "chicago", 1.8,
  "san-diego", "topeka", 1.4
)

costdf = costdf %>% mutate(cost = dtab * 90)

dict_from_df = function(df) {
  names(df) = c("x", "y", "z")
  dfpy = r_to_py(df)$set_index(c("x", "y"))
  dfdict = dfpy$to_dict()["z"]
  return(dfdict)
}

# this is a round about way to get a dictionary that has a tuple as key
# couldn't figure out how to write a list that would get converted to this dictionary
cost = dict_from_df(costdf %>% select(plants, mkts, cost))
cost

# create the model object
M = pyo$ConcreteModel()

# set the model parameters
M$plants = pyo$Set(initialize = c("seattle", "san-diego"))
M$mkts = pyo$Set(initialize = c("new-york", "chicago", "topeka"))

M$cap = pyo$Param(M$plants, initialize = list("seattle" = 350, "san-diego" = 600))
M$dem = pyo$Param(M$mkts, initialize = list("new-york" = 325, "chicago" = 300, "topeka" = 275))
M$cost = pyo$Param(M$plants, M$mkts, initialize = cost)

# define the model decision variables (shipment quantities between a plant and market)
M$x = pyo$Var(M$plants, M$mkts, bounds = tuple(0, NULL))

# supply from plants <= plant capacity
supply_rule = function(m, i) {
  supply = 0
  for (j in bi$list(m$mkts)) {
    supply = supply + m$x[tuple(i, j)]
  }
  expr = supply <= m$cap[[i]]
  return(expr)
}
M$supply_cons = pyo$Constraint(M$plants, rule = supply_rule)

# supply to a market meets its demand
demand_rule = function(m, j) {
  demand = 0
  for (i in bi$list(m$plants)) {
    demand = demand + m$x[tuple(i, j)]
  }
  expr = demand >= m$dem[[j]]
  return(expr)
}
M$demand_cons = pyo$Constraint(M$mkts, rule = demand_rule)  

# objective to minimize shipping cost
objective_rule = function(m) {
  totcost = 0
  for (i in bi$list(m$plants)) {
    for (j in bi$list(m$mkts)) {
      totcost = totcost + m$cost[tuple(i, j)] * m$x[tuple(i, j)]
    }
  }
  return(totcost)
}
M$objective = pyo$Objective(rule = objective_rule, sense = pyo$minimize)

# print model
cat(mdl_print("M"))

# solve using solver glpk
opt = pyo$SolverFactory('glpk')
results = opt$solve(M)

# get results into a dataframe
res_L = list()
k = 1
for(i in bi$list(M$plants)) {
    for(j in bi$list(M$mkts)) {
        res_L[[k]] = list(plant = i, mkt = j, qty = py_to_r(M$x[tuple(i, j)]()))
        k = k + 1
  }
}
res_df = bind_rows(res_L)
res_df


#'### Session Info
sessionInfo()
# pyomo version used is 5.7