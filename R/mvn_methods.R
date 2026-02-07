#' Generate Random Draws from a Multivariate Normal Mixture
#'
#' @description Draw random observations from a finite mixture of multivariate
#' normal components.
#' @export
#' @title rmvnmix
#' @name rmvnmix
#' @param n Number of observations.
#' @param alpha Length-\eqn{m} vector of component probabilities.
#' @param mu \eqn{d \times m} matrix whose \eqn{j}-th column is component mean.
#' @param sigma \eqn{d \times (d m)} block matrix of component covariance matrices.
#' @return An \eqn{n \times d} matrix of random draws.
rmvnmix <- function(n, alpha, mu, sigma) {
  m <- length(alpha)
  d <- nrow(mu)
  ind <- sample.int(m, size = n, replace = TRUE, prob = alpha)
  y <- matrix(0, nrow = n, ncol = d)

  for (j in seq_len(m)) {
    nj <- sum(ind == j)
    if (nj == 0L) {
      next
    }
    muj <- mu[, j]
    sigmaj <- sigma[, ((j - 1L) * d + 1L):(j * d), drop = FALSE]
    y[ind == j, ] <- mvtnorm::rmvnorm(nj, mean = muj, sigma = sigmaj)
  }

  y
}

#' Convert Half-Vectorized Covariance to Matrix
#'
#' @description Convert a lower-triangular half-vectorization
#' (vech order) into a symmetric covariance matrix.
#' @export
#' @title sigmavec2mat
#' @name sigmavec2mat
#' @param sigma.vec Numeric vector of length \eqn{d(d+1)/2}.
#' @param d Matrix dimension.
#' @return A \eqn{d \times d} symmetric matrix.
sigmavec2mat <- function(sigma.vec, d) {
  mvn_sigma_vech_to_mat(sigma_vec = sigma.vec, d = d)
}
