library(reticulate)

pd = import("pandas")

df = data.frame(x = c(1, 2, 3), y = c(10, 20, 30))

df2 = r_to_py(df)

df3 = df2$set_index("x")

py_to_r(df3)
