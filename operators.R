op = import("operator")

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

"<=.python.builtin.object" = function(a, b) {
  return(op$le(a, b))
}

">=.python.builtin.object" = function(a, b) {
  return(op$ge(a, b))
}

"==.python.builtin.object" = function(a, b) {
  return(op$eq(a, b))
}

