library(reticulate)

inspect = import("inspect")

fn1 = function(x, y) {
  z = x + y
  return(z)
}

fn1(3, 2)

fn1_py = r_to_py(fn1)

fn1_py(3, 2)

inspect$getargspec(fn1_py)
cat(inspect$getsource(fn1_py))

fn2_py = py_func(fn1)

inspect$getargspec(fn2_py)

