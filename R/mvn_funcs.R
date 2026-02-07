#' Generate Initial Values for PMLE
#' @description Generate initial values used by PMLE of multivariate normal
#' mixtures.
#' @export
#' @title mvnmixPMLEinit
#' @name mvnmixPMLEinit
#' @param y n by d matrix of data
#' @param ninits number of initial values to be generated
#' @param m The number of components in the mixture
#' @return A list with elements `alpha`, `mu`, and `sigma`.
mvnmixPMLEinit <- function(y, ninits = 1, m = 2) {
  y <- mvn_validate_y(y, allow_univariate = FALSE)
  d <- ncol(y)
  dsig <- d * (d + 1) / 2

  alpha <- matrix(stats::runif(m * ninits), nrow = m)
  alpha <- t(t(alpha) / colSums(alpha))

  mu <- matrix(0, nrow = d, ncol = m * ninits)
  variance <- matrix(0, nrow = d, ncol = m * ninits)

  for (i in seq_len(d)) {
    y0 <- y[, i]
    mu[i, ] <- stats::runif(m * ninits, min = min(y0), max = max(y0))
    variance_scale <- max(stats::var(y0), 1e-8)
    variance[i, ] <- stats::runif(m * ninits, min = 0.1, max = 2) * variance_scale
  }

  sigma <- mvn_draw_sigma_vech(variance = variance, corrmax = 0.4)

  list(
    alpha = alpha,
    mu = matrix(mu, nrow = d * m),
    sigma = matrix(sigma, nrow = dsig * m)
  )
}

#' Estimate PMLE for Multivariate Normal Mixtures
#' @description Estimates parameters of a finite mixture of multivariate normals
#' by penalized maximum log-likelhood functions.
#' @export
#' @title mvnmixPMLE
#' @name mvnmixPMLE
#' @param y n by d matrix of data
#' @param m The number of components in the mixture
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit The maximum number of iterations.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param binit The initial value of parameter vector that is included as a candidate parameter vector
#' @return A list containing estimated `coefficients`, `parlist`, `loglik`,
#'   `penloglik`, `aic`, `bic`, `postprobs`, `call`, and `m`.
#' @note \code{mvnmixPMLE} maximizes the penalized log-likelihood function
#' using the EM algorithm with combining short and long runs of EM steps as in Biernacki et al. (2003).
#' \code{mvnmixPMLE} first runs the EM algorithm from \code{ninits}\eqn{* 4m(1 + p)} initial values
#' with the convertence criterion \code{epsilon.short} and \code{maxit.short}.
#' Then, \code{mvnmixPMLE} uses \code{ninits} best initial values to run the EM algorithm
#' with the convertence criterion \code{epsilon} and \code{maxit}.
#' @references Alexandrovich, G. (2014)
#' A Note on the Article `Inference for Multivariate Normal Mixtures' by J. Chen and X. Tan
#' \emph{Journal of Multivariate Analysis}, \bold{129}, 245--248.
#'
#' Biernacki, C., Celeux, G. and Govaert, G. (2003)
#' Choosing Starting Values for the EM Algorithm for Getting the
#' Highest Likelihood in Multivariate Gaussian Mixture Models,
#' \emph{Computational Statistics and Data Analysis}, \bold{41}, 561--575.
#'
#' Boldea, O. and Magnus, J. R. (2009)
#' Maximum Likelihood Estimation of the Multivariate Normal Mixture Model,
#' \emph{Journal of the American Statistical Association},
#' \bold{104}, 1539--1549.
#'
#' Chen, J. and Tan, X. (2009)
#' Inference for Multivariate Normal Mixtures,
#' \emph{Journal of Multivariate Analysis}, \bold{100}, 1367--1383.
#'
#' McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley and Sons, Inc.
mvnmixPMLE <- function(y, m = 2,
                       ninits = 100, epsilon = 1e-08, maxit = 2000,
                       epsilon.short = 1e-02, maxit.short = 500, binit = NULL) {
  y <- mvn_validate_y(y, allow_univariate = FALSE)

  if (!is.numeric(m) || length(m) != 1L || m < 1L) {
    stop("m must be a positive integer.")
  }
  m <- as.integer(m)

  if (!is.numeric(ninits) || length(ninits) != 1L || ninits < 1L) {
    stop("ninits must be a positive integer.")
  }
  ninits <- as.integer(ninits)

  n <- nrow(y)
  d <- ncol(y)
  dsig <- d * (d + 1) / 2
  k_par <- mvn_param_count(m = m, d = d)

  var0 <- stats::var(y) * (n - 1) / n

  if (m == 1L) {
    mu <- colMeans(y)
    sigma <- var0[lower.tri(var0, diag = TRUE)]

    logdet_var0 <- as.numeric(determinant(var0, logarithm = TRUE)$modulus)
    loglik <- -(n / 2) * (d + d * log(2 * pi) + logdet_var0)
    penloglik <- loglik
    aic <- -2 * loglik + 2 * k_par
    bic <- -2 * loglik + log(n) * k_par

    parlist <- list(alpha = 1, mu = mu, sigma = sigma)
    coefficients <- c(alpha = 1, mu = mu, sigma = sigma)
    postprobs <- matrix(1, nrow = n, ncol = 1L)
  } else {
    ninits.short <- as.integer(max(1L, ninits * 10L * m * d))

    init <- mvnmixPMLEinit(y = y, ninits = ninits.short, m = m)
    b0 <- as.matrix(rbind(init$alpha, init$mu, init$sigma))

    if (!is.null(binit)) {
      if (length(binit) != nrow(b0)) {
        stop("binit has incompatible length.")
      }
      b0[, 1L] <- binit
    }

    an <- 1 / sqrt(n)
    var0vec <- var0[lower.tri(var0, diag = TRUE)]
    sigma0 <- rep(var0vec, m)
    mu0 <- numeric(m + 1L)

    out.short <- cppMVNmixPMLE(
      bs = b0, ys = y, mu0s = mu0, sigma0s = sigma0,
      m = m, an = an, maxit = maxit.short,
      ninits = ninits.short, tol = epsilon.short
    )

    components <- utils::head(order(out.short$penloglikset, decreasing = TRUE), ninits)
    b1 <- out.short$b[, components, drop = FALSE]

    out <- cppMVNmixPMLE(
      bs = b1, ys = y, mu0s = mu0, sigma0s = sigma0,
      m = m, an = an, maxit = maxit,
      ninits = ninits, tol = epsilon
    )

    index <- which.max(out$penloglikset)
    pars <- mvn_unpack_params(out$b, m = m, d = d, dsig = dsig, index = index)

    ordered <- mvn_order_components(
      alpha = pars$alpha,
      mu = pars$mu,
      sigma = pars$sigma,
      d = d,
      dsig = dsig
    )

    alpha <- ordered$alpha
    mu <- ordered$mu
    sigma <- ordered$sigma

    penloglik <- out$penloglikset[index]
    loglik <- out$loglikset[index]
    aic <- -2 * loglik + 2 * k_par
    bic <- -2 * loglik + log(n) * k_par

    postprobs <- matrix(out$post[, index], nrow = n)
    postprobs <- postprobs[, ordered$order, drop = FALSE]
    colnames(postprobs) <- paste0("comp.", seq_len(m))

    parlist <- list(alpha = alpha, mu = mu, sigma = sigma)
    coefficients <- unlist(parlist)
  }

  list(
    coefficients = coefficients,
    parlist = parlist,
    loglik = loglik,
    penloglik = penloglik,
    aic = aic,
    bic = bic,
    postprobs = postprobs,
    call = match.call(),
    m = m
  )
}

#' Compute the Max-Phi Objective for MEM Steps
#' @description Compute ordinary and penalized log-likelihood values from the
#' MEM algorithm at k = 1, 2, 3.
#' @title mvnmixMaxPhi
#' @name mvnmixMaxPhi
#' @param y n by d matrix of data
#' @param parlist The parameter estimates as a list containing alpha, mu, and sigma
#' in the form of (alpha = (alpha_1,...,alpha_m), mu = (mu_1',...,mu_m'),
#' sigma = (vech(sigma_1)',...,vech(sigma_m)')
#' @param an a term used for penalty function
#' @param tauset A set of initial tau value candidates
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @param verb Determines whether to print a message if an error occurs.
#' @param parallel Determines what percentage of available cores are used,
#'   represented by a number between 0 and 1.
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically create a new one for computation accordingly.
#' @return A list with `loglik`, `penloglik`, and `coefficient`.
mvnmixMaxPhi <- function(y, parlist, an, tauset = c(0.1, 0.3, 0.5),
                         ninits = 10, epsilon.short = 1e-02, epsilon = 1e-08,
                         maxit.short = 500, maxit = 2000,
                         verb = FALSE,
                         parallel = 0.75,
                         cl = NULL) {
  y <- mvn_validate_y(y, allow_univariate = TRUE)

  m <- length(parlist$alpha)
  d <- ncol(y)
  dsig <- d * (d + 1) / 2
  ninits.short <- as.integer(max(1L, ninits * 10L * m))

  loglik.all <- matrix(0, nrow = m * length(tauset), ncol = 3)
  penloglik.all <- matrix(0, nrow = m * length(tauset), ncol = 3)
  coefficient.all <- matrix(0, nrow = m * length(tauset), ncol = (1 + d + dsig) * (m + 1L))

  for (h in seq_len(m)) {
    for (t in seq_along(tauset)) {
      rowindex <- (t - 1L) * m + h
      tau <- tauset[t]
      result <- mvnmixMaxPhiStep(
        htaupair = c(h, tau), y = y, parlist = parlist, an = an,
        ninits = ninits, ninits.short = ninits.short,
        epsilon.short = epsilon.short, epsilon = epsilon,
        maxit.short = maxit.short, maxit = maxit,
        verb = verb
      )
      loglik.all[rowindex, ] <- result$loglik
      penloglik.all[rowindex, ] <- result$penloglik
      coefficient.all[rowindex, ] <- result$coefficient
    }
  }

  list(
    coefficient = as.vector(coefficient.all[which.max(loglik.all[, 3]), ]),
    loglik = apply(loglik.all, 2, max),
    penloglik = apply(penloglik.all, 2, max)
  )
}

#' Compute One (h, tau) Step for Max-Phi
#' @description Given a pair of h and tau and data, compute ordinary and
#' penalized log-likelihood values from MEM at k = 1, 2, 3.
#' @title mvnmixMaxPhiStep
#' @name mvnmixMaxPhiStep
#' @param htaupair A set of h and tau
#' @param y n by d matrix of data
#' @param parlist The parameter estimates as a list containing alpha, mu, and sigma
#' in the form of (alpha = (alpha_1,...,alpha_m), mu = (mu_1',...,mu_m'),
#' sigma = (vech(sigma_1)',...,vech(sigma_m)')
#' @param an a term used for penalty function
#' @param ninits The number of randomly drawn initial values.
#' @param ninits.short The number of candidates used to generate an initial phi, in short MEM
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @param verb Determines whether to print a message if an error occurs.
#' @return A list of phi, log-likelihood, and penalized log-likelihood resulting from MEM algorithm.
mvnmixMaxPhiStep <- function(htaupair, y, parlist, an,
                             ninits, ninits.short,
                             epsilon.short, epsilon,
                             maxit.short, maxit,
                             verb) {
  y <- mvn_validate_y(y, allow_univariate = TRUE)

  alpha0 <- parlist$alpha
  m <- length(alpha0)
  m1 <- m + 1L
  d <- ncol(y)
  dsig <- d * (d + 1) / 2

  h <- as.integer(htaupair[1])
  tau <- as.numeric(htaupair[2])

  if (h < 1L || h > m) {
    stop("h must be between 1 and m.")
  }

  mu0 <- parlist$mu
  mu0matrix <- matrix(mu0, nrow = d, ncol = m)
  mu0h <- c(-1e10, mu0matrix[1, ], 1e10)

  sigma0 <- parlist$sigma
  sigma0h <- c(sigma0[1:(h * dsig)], sigma0[((h - 1L) * dsig + 1L):(m * dsig)])

  init <- mvnmixPhiInit(y = y, parlist = parlist, h = h, tau = tau, ninits = ninits.short)
  b0 <- as.matrix(rbind(init$alpha, init$mu, init$sigma))

  out.short <- cppMVNmixPMLE(
    bs = b0, ys = y, mu0s = mu0h, sigma0s = sigma0h,
    m = m1, an = an, maxit = maxit.short,
    ninits = ninits.short, tol = epsilon.short,
    tau = tau, h = h, k = 1L
  )

  components <- utils::head(order(out.short$penloglikset, decreasing = TRUE), ninits)
  if (verb && any(out.short$notcg)) {
    cat(sprintf("non-convergence rate at short-EM = %.3f\n", mean(out.short$notcg)))
  }

  out <- cppMVNmixPMLE(
    bs = out.short$b[, components, drop = FALSE], ys = y,
    mu0s = mu0h, sigma0s = sigma0h,
    m = m1, an = an, maxit = maxit,
    ninits = ninits, tol = epsilon,
    tau = tau, h = h, k = 1L
  )

  index <- which.max(out$penloglikset)
  pars <- mvn_unpack_params(out$b, m = m1, d = d, dsig = dsig, index = index)
  ordered <- mvn_order_components(pars$alpha, pars$mu, pars$sigma, d = d, dsig = dsig)

  sigma0h.matrix <- matrix(sigma0h, nrow = dsig, ncol = m1)
  sigma0h <- c(sigma0h.matrix[, ordered$order, drop = FALSE])

  b <- as.matrix(c(ordered$alpha, ordered$mu, ordered$sigma))

  loglik <- numeric(3)
  penloglik <- numeric(3)
  loglik[1] <- out$loglikset[[index]]
  penloglik[1] <- out$penloglikset[[index]]

  alpha <- ordered$alpha
  mu <- ordered$mu
  sigma <- ordered$sigma

  for (step_k in 2:3) {
    out <- cppMVNmixPMLE(
      bs = b, ys = y, mu0s = mu0h, sigma0s = sigma0h,
      m = m1, an = an, maxit = 2L,
      ninits = 1L, tol = epsilon,
      tau = tau, h = h, k = step_k
    )

    b <- out$b
    pars <- mvn_unpack_params(b, m = m1, d = d, dsig = dsig, index = 1L)

    alpha <- pars$alpha
    mu <- pars$mu
    sigma <- pars$sigma

    loglik[step_k] <- out$loglikset[[1]]
    penloglik[step_k] <- out$penloglikset[[1]]

    if (mvn_is_singular_state(alpha = alpha, sigma = sigma, d = d, tol = 1e-6)) {
      loglik[step_k] <- -Inf
      penloglik[step_k] <- -Inf
      break
    }

    ordered <- mvn_order_components(alpha, mu, sigma, d = d, dsig = dsig)
    alpha <- ordered$alpha
    mu <- ordered$mu
    sigma <- ordered$sigma
    sigma0h <- c(matrix(sigma0h, nrow = dsig, ncol = m1)[, ordered$order, drop = FALSE])

    b <- as.matrix(c(alpha, mu, sigma))
  }

  list(
    coefficient = as.matrix(c(alpha, mu, sigma)),
    loglik = loglik,
    penloglik = penloglik
  )
}

#' Generate Initial Phi Candidates for MEM
#' @description Generates parameter candidates used by the modified EM test for
#' mixtures of multivariate normals.
#' @title mvnmixPhiInit
#' @name mvnmixPhiInit
#' @param y n by d matrix of data
#' @param parlist The parameter estimates as a list containing alpha, mu, and sigma
#' in the form of (alpha = (alpha_1,...,alpha_m),
#' mu = (mu_1',...,mu_m'), sigma = (vech(sigma_1)',...,vech(sigma_m)')
#' @param h h used as index for pivoting
#' @param tau Tau used to split the h-th component
#' @param ninits number of initial values to be generated
#' @return A list with elements `alpha`, `mu`, and `sigma`.
mvnmixPhiInit <- function(y, parlist, h, tau, ninits = 1) {
  y <- mvn_validate_y(y, allow_univariate = TRUE)

  d <- ncol(y)
  dsig <- d * (d + 1) / 2
  mu0 <- parlist$mu
  sigma0 <- parlist$sigma
  alpha0 <- parlist$alpha
  m <- length(alpha0)

  mu0matrix <- matrix(mu0, nrow = d, ncol = m)
  sigma0matrix <- matrix(sigma0, nrow = dsig, ncol = m)
  y1 <- y[, 1]
  mu01 <- mu0matrix[1, ]

  if (m >= 2L) {
    mid <- (mu01[1:(m - 1L)] + mu01[2:m]) / 2
    lb0 <- c(min(y1), mid)
    ub0 <- c(mid, max(y1))
    lb <- c(lb0[1:h], lb0[h:m])
    ub <- c(ub0[1:h], ub0[h:m])
  } else {
    lb <- c(min(y1), min(y1))
    ub <- c(max(y1), max(y1))
  }

  mu <- matrix(0, nrow = d, ncol = (m + 1L) * ninits)
  variance <- matrix(0, nrow = d, ncol = (m + 1L) * ninits)

  ninits1 <- floor(ninits / 2)
  ninits2 <- ninits - ninits1

  mu[1, ] <- stats::runif((m + 1L) * ninits, min = lb, max = ub)

  sigma01 <- sigma0matrix[1, ]
  sigma.1.hyp <- c(sigma01[1:h], sigma01[h:m])
  sigma.1.hyp <- pmax(sigma.1.hyp, 1e-8)
  variance[1, ] <- stats::runif((m + 1L) * ninits, min = sigma.1.hyp * 0.25, max = sigma.1.hyp * 2)

  for (i in 2:d) {
    y.i <- y[, i]
    mu0i <- mu0matrix[i, ]
    sigma0i <- sigma0matrix[(i - 1L) * d - (i - 1L) * (i - 2L) / 2 + 1L, ]

    mu.i.hyp <- c(mu0i[1:h], mu0i[h:m])
    sigma.i.hyp <- pmax(c(sigma0i[1:h], sigma0i[h:m]), 1e-8)
    sd.i.hyp <- sqrt(sigma.i.hyp)

    mu.i1 <- stats::runif((m + 1L) * ninits1, min = mu.i.hyp - sd.i.hyp, max = mu.i.hyp + sd.i.hyp)
    mu.i2 <- stats::runif((m + 1L) * ninits2, min = min(y.i), max = max(y.i))

    mu[i, ] <- c(mu.i1, mu.i2)
    variance[i, ] <- stats::runif((m + 1L) * ninits, min = sigma.i.hyp * 0.25, max = sigma.i.hyp * 2)
  }

  sigma <- mvn_draw_sigma_vech(variance = variance, corrmax = 0.4)

  alpha.hyp <- c(alpha0[1:h], alpha0[h:m])
  alpha.hyp[h:(h + 1L)] <- c(alpha.hyp[h] * tau, alpha.hyp[h + 1L] * (1 - tau))
  alpha <- matrix(rep.int(alpha.hyp, ninits), nrow = m + 1L)

  list(
    alpha = alpha,
    mu = matrix(mu, nrow = d * (m + 1L)),
    sigma = matrix(sigma, nrow = dsig * (m + 1L))
  )
}
