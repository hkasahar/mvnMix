mvnmixLRT <- function(y, m = 2, an = 1, tauset = c(0.1, 0.3, 0.5),
                      ninits = 10,
                      crit.method = c("asy", "boot", "none"), nbtsp = 199,
                      cl = NULL,
                      parallel = 0.75,
                      LRT.penalized = FALSE) {
  y <- mvn_validate_y(y, allow_univariate = FALSE)
  crit.method <- match.arg(crit.method)

  pmle.result <- mvnmixPMLE(y = y, m = m, ninits = ninits)
  loglik0 <- pmle.result$loglik

  pmle1.result <- mvnmixPMLE(y = y, m = m + 1L, ninits = ninits)
  loglik1 <- pmle1.result$loglik

  lrt <- 2 * (loglik1 - loglik0)
  if (LRT.penalized) {
    lrt <- 2 * (loglik1 - loglik0)
  }

  a <- list(LRT = lrt, ll0 = loglik0, ll1 = loglik1,
            call = match.call(), m = m, label = "MEMtest")

  class(a) <- "normalregMix"
  a
}

#' Perform Modified EM Test for Mixture Order
#' @export
#' @title mvnmixMEMtest
#' @name mvnmixMEMtest
#' @param y n by d matrix of data
#' @param m The number of components in the mixture defined by a null hypothesis, m_0
#' @param tauset A set of initial tau value candidates
#' @param an a term used for penalty function
#' @param ninits The number of randomly drawn initial values.
#' @param crit.method Method used to compute critical values: one of
#'   \code{"none"}, \code{"asy"}, and \code{"boot"}.
#' @param nbtsp The number of bootstrap observations.
#' @param cl Cluster used for parallelization.
#' @param parallel Share of available cores used for bootstrap simulation.
#' @param LRT.penalized Determines whether penalized likelihood is used in the
#'   alternative model objective.
#' @param nrep Number of simulation draws used in asymptotic critical values.
#' @param ninits.crit Number of random starts used in cone-projection steps.
#' @return A list containing `emstat`, `pvals`, `crit`, fitted null-model
#'   outputs, and metadata.
#' @examples
#' data(faithful)
#' y <- as.matrix(faithful[, c("eruptions", "waiting")])
#' mvnmixMEMtest(y = y, m = 1, crit.method = "none", ninits = 3)
mvnmixMEMtest <- function(y, m = 2, an = 1, tauset = c(0.1, 0.3, 0.5),
                          ninits = 10,
                          crit.method = c("asy", "boot", "none"), nbtsp = 199,
                          cl = NULL,
                          parallel = 0.75,
                          LRT.penalized = FALSE,
                          nrep = 10000, ninits.crit = 25) {
  y <- mvn_validate_y(y, allow_univariate = FALSE)
  n <- nrow(y)
  d <- ncol(y)
  crit.method <- match.arg(crit.method)

  pmle.result <- mvnmixPMLE(y = y, m = m, ninits = ninits)
  loglik0 <- pmle.result$loglik

  sigma.matrix <- matrix(pmle.result$parlist$sigma, ncol = m)
  diag_idx <- mvn_diag_vech_index(d)
  sigma.diags <- sigma.matrix[diag_idx, , drop = FALSE]
  min_diag <- min(sigma.diags)

  if (!is.finite(min_diag) || min_diag <= 0) {
    sc <- 1
  } else {
    sc <- sqrt(0.1 / min_diag)
  }

  y_scaled <- y * sc
  parlist0 <- pmle.result$parlist
  parlist0$mu <- parlist0$mu * sc
  parlist0$sigma <- parlist0$sigma * sc * sc

  par1 <- mvnmixMaxPhi(
    y = y_scaled, parlist = parlist0,
    an = an, tauset = tauset, ninits = ninits,
    parallel = 0, cl = cl
  )

  emstat <- 2 * (par1$loglik - loglik0 + log(sc) * n * d)
  if (LRT.penalized) {
    emstat <- 2 * (par1$penloglik - loglik0 + log(sc) * n * d)
  }

  if (crit.method == "asy") {
    result <- mvnmixCrit(y = y_scaled, parlist = parlist0, values = emstat,
                         nrep = nrep, ninits.crit = ninits.crit)
  } else if (crit.method == "boot") {
    result <- mvnmixCritBoot(
      y = y_scaled, an = an, parlist = parlist0, values = emstat,
      ninits = ninits, nbtsp = nbtsp, parallel = parallel, cl = cl,
      LRT.penalized = LRT.penalized
    )
  } else {
    result <- list(crit = rep(NA_real_, 3), pvals = rep(NA_real_, 3))
  }

  a <- list(
    emstat = emstat, pvals = result$pvals, crit = result$crit, crit.method = crit.method,
    parlist = pmle.result$parlist, ll0 = loglik0, ll1 = par1$loglik,
    aic = pmle.result$aic, bic = pmle.result$bic, postprobs = pmle.result$postprobs,
    call = match.call(), m = m, label = "MEMtest"
  )

  class(a) <- "normalregMix"
  a
}

#' Compute Bootstrap Critical Values for MEM Test
#' @description Computes bootstrap critical values for the modified EM test.
#' @export
#' @title mvnmixCritBoot
#' @name mvnmixCritBoot
#' @param y n by d matrix of data
#' @param an Penalty-weight parameter used in PMLE updates.
#' @param parlist Parameter list with `alpha`, `mu`, and `sigma`.
#' @param values Optional vector of observed statistics for p-value calculation.
#' @param ninits The number of initial candidates to be generated.
#' @param nbtsp The number of bootstrap observations.
#' @param parallel Determines what percentage of available cores are used.
#' @param cl Cluster used for parallelization (optional).
#' @param LRT.penalized Whether penalized log-likelihood is used in bootstrap
#'   test-statistic updates.
#' @return A list with elements `crit` and `pvals`.
mvnmixCritBoot <- function(y, an = 1, parlist, values = NULL, ninits = 10,
                           nbtsp = 199, parallel = 0.75, cl = NULL,
                           LRT.penalized = FALSE) {
  y <- mvn_validate_y(y, allow_univariate = FALSE)

  n <- nrow(y)
  d <- ncol(y)
  alpha <- parlist$alpha
  mu <- parlist$mu
  sigma <- parlist$sigma
  m <- length(alpha)

  mu.mat <- matrix(mu, nrow = d, ncol = m)
  sigma.mat <- mvn_sigma_vec_to_blocks(sigma = sigma, d = d, m = m)

  if (m == 1L) {
    ybset <- mvtnorm::rmvnorm(nbtsp * n, mean = mu, sigma = sigma.mat)
  } else {
    ybset <- rmvnmix(nbtsp * n, alpha = alpha, mu = mu.mat, sigma = sigma.mat)
  }
  ybset <- array(ybset, dim = c(n, nbtsp, d))

  out <- NULL
  num.cores <- max(1L, floor(parallel::detectCores() * parallel))

  if (num.cores > 1L) {
    if (is.null(cl)) {
      cl <- parallel::makeCluster(num.cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
    }

    parallel::clusterSetRNGStream(cl = cl, iseed = 123456)
    out <- parallel::parLapply(cl, seq_len(nbtsp), function(j) {
      mvnmixMEMtest(
        y = ybset[, j, ], m = m, parallel = 0,
        an = an, ninits = ninits, crit.method = "none",
        LRT.penalized = LRT.penalized
      )
    })
  } else {
    out <- lapply(seq_len(nbtsp), function(j) {
      mvnmixMEMtest(
        y = ybset[, j, ], m = m, parallel = 0,
        an = an, ninits = ninits, crit.method = "none",
        LRT.penalized = LRT.penalized
      )
    })
  }

  emstat.b <- sapply(out, "[[", "emstat")
  emstat.b <- t(apply(emstat.b, 1, sort))

  q <- ceiling(nbtsp * c(0.90, 0.95, 0.99))
  crit <- emstat.b[, q, drop = FALSE]

  pvals <- NULL
  if (!is.null(values)) {
    pvals <- rowMeans(emstat.b > values)
  }

  list(crit = crit, pvals = pvals)
}
