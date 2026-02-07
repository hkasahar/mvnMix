test_that("mvnmixPMLE rejects univariate input", {
  y <- matrix(stats::rnorm(60), ncol = 1)
  expect_error(mvnmixPMLE(y, m = 1), "more than one columns")
})

test_that("mvnmixPMLE runs in higher dimension with finite outputs", {
  set.seed(202601)
  n <- 120
  d <- 6

  y <- mvtnorm::rmvnorm(n = n, mean = rep(0, d), sigma = diag(d))

  fit <- mvnmixPMLE(
    y,
    m = 2,
    ninits = 3,
    epsilon = 1e-6,
    maxit = 300,
    epsilon.short = 1e-3,
    maxit.short = 80
  )

  expect_true(is.finite(fit$loglik))
  expect_true(is.finite(fit$penloglik))
  expect_equal(dim(fit$postprobs), c(n, 2))
  expect_equal(sum(fit$parlist$alpha), 1, tolerance = 1e-6)
})
