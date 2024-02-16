test_that("coxtp fits models correctly", {
  data(ExampleData)
  z <- ExampleData$z
  time <- ExampleData$time
  event <- ExampleData$event
  lambda <- c(0, 1)
  
  fit <- coxtp(event = event, z = z, time = time, lambda = lambda)
  
  expect_true("coxtp" %in% class(fit), "fit should have class 'coxtp'")
  
  expect_length(fit, length(lambda) + 1)
})