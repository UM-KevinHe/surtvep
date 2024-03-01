library(surtvep)
library(testthat)
test_that("The location of the self-specified knots", code = {
  data(ExampleData)
  z     <- ExampleData$z
  time  <- ExampleData$time
  event <- ExampleData$event
  fit.tp <- coxtp(event = event, z = z, time = time, lambda = 100, penalty = "Smooth-spline", nspline = 7, knots = c(0.5, 1.0, 2.0))
  expect_identical(fit.tp$lambda1$internal.knots, c(0.5, 1, 2))
})
