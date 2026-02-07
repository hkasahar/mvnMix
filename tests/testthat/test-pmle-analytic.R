test_that("mvnmixPMLE m=1 matches analytical moments and loglik", {
  y <- rbind(
    c(0.5, -0.4),
    c(1.2, 0.3),
    c(0.1, 1.1),
    c(-0.8, 0.7),
    c(0.9, -0.2)
  )

  fit <- mvnmixPMLE(y, m = 1)

  n <- nrow(y)
  d <- ncol(y)
  sigma0 <- stats::cov(y) * (n - 1) / n
  expected_loglik <- -(n / 2) * (
    d + d * log(2 * pi) +
      as.numeric(determinant(sigma0, logarithm = TRUE)$modulus)
  )

  expect_equal(fit$parlist$alpha, 1)
  expect_equal(unname(fit$parlist$mu), unname(colMeans(y)), tolerance = 1e-10)
  expect_equal(
    unname(fit$parlist$sigma),
    unname(sigma0[lower.tri(sigma0, diag = TRUE)]),
    tolerance = 1e-10
  )
  expect_equal(fit$loglik, expected_loglik, tolerance = 1e-10)
})
