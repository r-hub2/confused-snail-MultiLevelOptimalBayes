## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(MultiLevelOptimalBayes)

## ----eval = FALSE-------------------------------------------------------------
# mlob(
#   formula,
#   data,
#   group,
#   balancing.limit = 0.2,
#   conf.level = 0.95,
#   jackknife = FALSE,
#   punish.coeff = 2,
#   ...
# )

## ----eval= FALSE--------------------------------------------------------------
# result_iris <- mlob(
#   Sepal.Length ~ Sepal.Width + Petal.Length,
#   data = iris,
#   group = "Species",
#   conf.level = 0.99,
#   jackknife = FALSE
# )
# 
# summary(result_iris)

## ----eval= FALSE--------------------------------------------------------------
# result_chick <- mlob(
#   weight ~ Diet,
#   data = ChickWeight,
#   group = "Time",
#   punish.coeff = 1.5,
#   jackknife = FALSE
# )
# 
# print(result_chick)
# summary(result_chick)

## ----eval= FALSE--------------------------------------------------------------
# result_mtcars <- mlob(
#   mpg ~ hp + wt + am + hp:wt + hp:am,
#   data = mtcars,
#   group = "cyl",
#   balancing.limit = 0.35
# )
# 
# summary(result_mtcars)

