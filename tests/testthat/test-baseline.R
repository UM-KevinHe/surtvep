test_that("baseline function computes correctly", {
  data(ExampleData)
  z <- ExampleData$z
  time <- ExampleData$time
  event <- ExampleData$event
  
  fit <- coxtv(event = event, z = z, time = time)
  base.est <- baseline(fit)
  
  # Check if 'base.est' has the expected class, for example, a data.frame
  expect_true("baseline" == class(base.est), "base.est should be a baseline class")
  
  # Check if 'base.est' contains expected columns, for example: 'time' and 'baseline'
  expect_true(all(c("time", "hazard", "cumulHaz") %in% names(base.est)),
              "base.est should contain 'time', 'hazard' and 'cumulHaz' columns")
})
