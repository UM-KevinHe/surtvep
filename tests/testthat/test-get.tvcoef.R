test_that("get.tvcoef returns correct structure and dimensions", {
  data(ExampleData)
  z <- ExampleData$z
  time <- ExampleData$time
  event <- ExampleData$event
  
  fit <- coxtv(event = event, z = z, time = time, degree = 2)
  
  coef <- get.tvcoef(fit)
  
  expect_true(is.matrix(coef), info = "coef should be a matrix")
  
  expect_equal(dim(coef), c(length(unique(time)), ncol(z)), 
               info = "Dimensions of coef should match length of unique times by number of covariates")
  
  expect_false(any(is.na(coef)), info = "Coefficient matrix should not contain NA values")
})
