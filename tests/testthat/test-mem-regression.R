test_that("mvnmixMEMtest regression: emstat stays on reference path", {
  ref <- readRDS(testthat::test_path("testdata", "memstat_reference.rds"))

  set.seed(ref$seed)
  y <- rmvnmix(
    n = ref$n,
    alpha = ref$alpha,
    mu = ref$mu,
    sigma = ref$sigma
  )

  out <- mvnmixMEMtest(
    y = y,
    m = ref$m_null,
    an = ref$an,
    tauset = ref$tauset,
    ninits = ref$ninits,
    crit.method = "none"
  )

  expect_equal(as.numeric(out[["emstat"]]), ref$emstat, tolerance = 1e-6)
})
