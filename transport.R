library(reticulate)
library(dplyr)

bi = import_builtins()

pyo = import("pyomo.environ", convert = FALSE)
source("operators.R")

plants = c("seattle", "san-diego")
mkts = c("new-york", "chicago", "topeka")

cap = list("seattle" = 350, "san-diego" = 600)
dem = list("new-york" = 325, "chicago" = 300, "topeka" = 275)

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

M = pyo$ConcreteModel()
M$plants = pyo$Set(initialize = c("seattle", "san-diego"))
M$mkts = pyo$Set(initialize = c("new-york", "chicago", "topeka"))

M$cap = pyo$Param(M$plants, initialize = list("seattle" = 350, "san-diego" = 600))
M$dem = pyo$Param(M$mkts, initialize = list("new-york" = 325, "chicago" = 300, "topeka" = 275))
M$cost = pyo$Param(M$plants, M$mkts, initialize = cost)

M$x = pyo$Var(M$plants, M$mkts, bounds = tuple(0, NULL))

supply_rule = function(m, i) {
  supply = 0
  for (j in bi$list(m$mkts)) {
#    print(j)
    supply = supply + m$x[tuple(i, j)]
  }
  expr = supply <= m$cap[[i]]
  return(expr)
}
M$supply_cons = pyo$Constraint(M$plants, rule = supply_rule)

demand_rule = function(m, j) {
  demand = 0
  for (i in bi$list(m$plants)) {
    demand = demand + m$x[tuple(i, j)]
  }
  expr = demand >= m$dem[[j]]
  return(expr)
}
M$demand_cons = pyo$Constraint(M$mkts, rule = demand_rule)  

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

opt = pyo$SolverFactory('glpk')
results = opt$solve(M)

M$x$display()

for(j in bi$list(M$mkts)) {
  print(j)
}
