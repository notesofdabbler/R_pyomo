#' ---
#' title: "Using Pyomo from R through the magic of Reticulate"
#' output: html_document
#' author: Shankar Vaidyaraman
#' ---
#'
#' [Pyomo](http://www.pyomo.org/) is a python based open-source package for modeling optimization
#' problems. It makes it easy to represent optimization problems and can send it to different solvers
#' (both open-source and commercial) to solve the problem and return the results in python. The advantage of
#' pyomo compared to commercial software such as [GAMS](https://www.gams.com/) and [AMPL](https://ampl.com/) is the ability to code using standard python
#' syntax (with some modifications for pyomo constructs). Another open source package for modeling optimization problems is
#' [JuMP](https://jump.dev/JuMP.jl/v0.19.0/index.html) in Julia language. 
#' 
#' My goal in this blog is to see how far I can get in terms of using Pyomo from R using the 
#' [reticulate](https://rstudio.github.io/reticulate/) package. The simplest option would be to 
#' develop the model in pyomo and call it from R using reticulate. While that is a simple option
#' and will work, it still requires writing the pyomo model in python. I want to use reticulate to be
#' able to write the pyomo model using R. This blog catalogs my learnings by trying to develop a few
#' optimization models in pyomo from R.
#' 
#' 
#'## Set-up
#+ set-up
# specify the python to use
reticulate::use_python("/Users/shanki/anaconda3/bin/python")
# load libraries
library(reticulate)
library(dplyr)
library(tidyr)
library(ggplot2)
# Import pyomo
pyo = import("pyomo.environ", convert = FALSE)
#' When trying different examples, I saw some issues when not using the `convert = FALSE` option and so started using
#' this option when importing pyomo
#' 
#' ## Rosenbrock Problem (Unconstrained Nonlinear Optimization)
#' The first example is a unconstrained nonlinear minimization of 
#' [rosenbrock function](https://en.wikipedia.org/wiki/Rosenbrock_function)
#'$$
#' \min \;\; (x - 1)^2 + (y - x^2)^2
#' $$ 
#' 
#' Since arithmetic operators in pyomo are overloaded to generate pyomo expressions, I wrote
#' functions for [arithmetic operators] based on this [discussion thread]() and sourced them here.
#' 
#+
source("operators.R")
#' 