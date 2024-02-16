test_that("cv.coxtp returns correct structure and valid results", {
  data(ExampleData)
  z <- ExampleData$z
  time <- ExampleData$time
  event <- ExampleData$event
  lambda <- c(0.1, 1)
  
  fit <- cv.coxtp(event = event, z = z, time = time, lambda = lambda, nfolds = 5)
  
  expect_true(inherits(fit, "cv.coxtp"), "fit should inherit from 'cv.coxtp'")
  
  expect_true(is.list(fit), "fit should be a list")
  expect_true(all(c("model.cv", "lambda", "cve", "lambda.min") %in% names(fit)),
              "fit should contain 'model.cv', 'lambda', 'cve', 'lambda.min'")
  
  expect_true(inherits(fit$model.cv, "coxtp"), "'model.cv' should be a 'coxtp' object")
  
  expect_equal(fit$lambda, lambda, info = "fit$lambda should match the input lambda values")
  
  expect_true(is.numeric(fit$cve) && length(fit$cve) == length(lambda),
              "'cve' should be a numeric vector with length equal to 'lambda'")
  
  expect_true(is.numeric(fit$lambda.min) && length(fit$lambda.min) == 1,
              "'lambda.min' should be a single numeric value")
  expect_true(fit$lambda.min %in% lambda, "'lambda.min' should be one of the input lambda values")
  
  expect_true(all(fit$cve >= 0), "'cve' should contain non-negative values")
})
