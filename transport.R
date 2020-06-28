library(reticulate)
library(dplyr)

source("operators.R")

pyo = import("pyomo.environ")

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

model = pyo$ConcreteModel()
model$x = pyo$Var(plants, mkts, bounds = tuple(0, NULL))

supply_rule = function(m, i) {
  supply = 0
  for (j in mkts) {
    supply = supply + m$x[tuple(i, j)]
  }
  expr = supply <= cap[[i]]
  return(expr)
}
model$supply_cons = pyo$Constraint(plants, rule = supply_rule)

demand_rule = function(m, j) {
  demand = 0
  for (i in plants) {
    demand = demand + m$x[tuple(i, j)]
  }
  expr = demand >= dem[[j]]
  return(expr)
}
model$demand_cons = pyo$Constraint(mkts, rule = demand_rule)  

objective_rule = function(m) {
  totcost = 0
  for (i in plants) {
    for (j in mkts) {
      totcost = totcost + cost[tuple(i, j)] * model$x[tuple(i, j)]
    }
  }
  return(totcost)
}
model$objective = pyo$Objective(rule = objective_rule, sense = pyo$minimize)

opt = pyo$SolverFactory('glpk')
results = opt$solve(model)

model$x$display()
