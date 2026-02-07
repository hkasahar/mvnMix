# Internal utilities for parameter reshaping, covariance handling, and
# numerically stable linear algebra used across mvnMix estimators/tests.

mvn_validate_y <- function(y, allow_univariate = TRUE) {
  y <- as.matrix(y)
  if (!is.numeric(y)) {
    stop("y must be a numeric matrix.")
  }
  if (nrow(y) < 2L) {
    stop("y must have at least two rows.")
  }
  if (!allow_univariate && ncol(y) < 2L) {
    stop("y must have more than one columns.")
  }
  if (any(!is.finite(y))) {
    stop("y must contain only finite values.")
  }
  y
}

mvn_diag_vech_index <- function(d) {
  cumsum(c(1L, seq.int(d, 2L)))
}

mvn_param_count <- function(m, d) {
  dsig <- d * (d + 1L) / 2L
  m - 1L + m * d + m * dsig
}

mvn_unpack_params <- function(b, m, d, dsig, index = 1L) {
  list(
    alpha = b[seq_len(m), index],
    mu = b[(m + 1L):(m + m * d), index],
    sigma = b[(m + m * d + 1L):(m + m * d + m * dsig), index]
  )
}

mvn_order_components <- function(alpha, mu, sigma, d, dsig) {
  mu_mat <- matrix(mu, nrow = d, ncol = length(alpha))
  sigma_mat <- matrix(sigma, nrow = dsig, ncol = length(alpha))
  ord <- order(mu_mat[1L, ])
  list(
    alpha = alpha[ord],
    mu = c(mu_mat[, ord, drop = FALSE]),
    sigma = c(sigma_mat[, ord, drop = FALSE]),
    order = ord
  )
}

mvn_sigma_vech_to_mat <- function(sigma_vec, d) {
  sigma <- diag(d)
  sigma[lower.tri(sigma, diag = TRUE)] <- sigma_vec
  sigma <- sigma + t(sigma)
  diag(sigma) <- diag(sigma) / 2
  sigma
}

mvn_sigma_vec_to_blocks <- function(sigma, d, m) {
  dsig <- d * (d + 1L) / 2L
  out <- matrix(0, nrow = d, ncol = d * m)
  for (j in seq_len(m)) {
    idx <- ((j - 1L) * dsig + 1L):(j * dsig)
    out[, ((j - 1L) * d + 1L):(j * d)] <- mvn_sigma_vech_to_mat(sigma[idx], d)
  }
  out
}

mvn_is_singular_state <- function(alpha, sigma, d, tol = 1e-6) {
  m <- length(alpha)
  dsig <- d * (d + 1L) / 2L
  if (any(!is.finite(alpha)) || any(alpha < tol)) {
    return(TRUE)
  }
  detsigma <- numeric(m)
  for (j in seq_len(m)) {
    idx <- ((j - 1L) * dsig + 1L):(j * dsig)
    sigma_j <- mvn_sigma_vech_to_mat(sigma[idx], d)
    detsigma[j] <- determinant(sigma_j, logarithm = TRUE)$modulus
  }
  any(!is.finite(detsigma)) || any(detsigma < log(tol))
}

mvn_draw_sigma_vech <- function(variance, corrmax = 0.4) {
  d <- nrow(variance)
  width <- ncol(variance)
  dsig <- d * (d + 1L) / 2L
  sigma <- matrix(0, nrow = dsig, ncol = width)

  sigma[1L, ] <- variance[1L, ]
  if (d >= 2L) {
    sigma[2:d, ] <-
      sqrt(t(t(variance[2:d, , drop = FALSE]) * variance[1L, ])) *
      matrix(stats::runif((d - 1L) * width, min = -corrmax, max = corrmax), nrow = d - 1L)
  }

  if (d >= 3L) {
    for (i in seq_len(d - 2L)) {
      row_diag <- i * d - i * (i - 1L) / 2L + 1L
      sigma[row_diag, ] <- variance[i + 1L, ]
      row_off <- (row_diag + 1L):((i + 1L) * d - i * (i + 1L) / 2L)
      sigma[row_off, ] <-
        sqrt(t(t(variance[(i + 2L):d, , drop = FALSE]) * variance[i + 1L, ])) *
        matrix(
          stats::runif((d - i - 1L) * width, min = -corrmax, max = corrmax),
          nrow = d - i - 1L
        )
    }
  }

  sigma[dsig, ] <- variance[d, ]
  sigma
}

mvn_chol_with_jitter <- function(mat, jitter = 1e-10, max_tries = 8L) {
  m <- (mat + t(mat)) / 2
  d <- nrow(m)
  diag_add <- 0
  for (k in 0:max_tries) {
    attempt <- tryCatch(chol(m + diag(diag_add, d)), error = function(e) NULL)
    if (!is.null(attempt)) {
      return(list(chol = attempt, mat = m + diag(diag_add, d), jitter = diag_add))
    }
    diag_add <- jitter * (10^k)
  }
  stop("Matrix is numerically singular even after diagonal regularization.")
}

mvn_stable_solve <- function(A, B = NULL, jitter = 1e-10, max_tries = 8L) {
  A <- (A + t(A)) / 2
  d <- nrow(A)
  diag_add <- 0
  for (k in 0:max_tries) {
    A_reg <- A + diag(diag_add, d)
    out <- tryCatch(
      {
        if (is.null(B)) {
          solve(A_reg)
        } else {
          solve(A_reg, B)
        }
      },
      error = function(e) NULL
    )
    if (!is.null(out)) {
      return(out)
    }
    diag_add <- jitter * (10^k)
  }
  stop("Unable to solve linear system after regularization.")
}

mvn_regularize_psd <- function(mat, eps = 1e-12) {
  mat <- (mat + t(mat)) / 2
  eig <- eigen(mat, symmetric = TRUE)
  vals <- pmax(eig$values, eps * max(1, abs(eig$values[1L])))
  eig$vectors %*% (vals * t(eig$vectors))
}
