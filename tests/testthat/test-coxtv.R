test_that("coxtv function works correctly", {
  # Assuming ExampleData is available in your package or test environment
  data(ExampleData)
  z <- ExampleData$z
  time <- ExampleData$time
  event <- ExampleData$event
  fit <- coxtv(event = event, z = z, time = time)
  
  expect_true(is.list(fit)) # Check if fit is a list
  expect_s3_class(fit, "coxtv") # If coxtv is expected to return a coxph object
})