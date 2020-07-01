library(reticulate)
library(dplyr)

bi = import_builtins()

pyo = import("pyomo.environ", convert = FALSE)
source("operators.R")

# plants
plants = c("seattle", "san-diego")
# markets
mkts = c("new-york", "chicago", "topeka")

# capacity of plants
cap = list("seattle" = 350, "san-diego" = 600)
# demand of markets
dem = list("new-york" = 325, "chicago" = 300, "topeka" = 275)

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

cost = dict_from_df(costdf %>% select(plants, mkts, cost))
dict_from_df(costdf)

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

# solve using solver glpk
opt = pyo$SolverFactory('glpk')
results = opt$solve(M)

M$x$display()

res_L = list()
k = 0
for(i in bi$list(M$plants)) {
    for(j in bi$list(M$mkts)) {
        res_L[[k]] = list(plant = i, mkt = j, qty = py_to_r(M$x[tupe(i, j)]()))
  }
}
res_df = bind_rows(res_L)
res_df
