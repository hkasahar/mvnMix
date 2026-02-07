#' Asymptotic critical values for multivariate normal mixture MEM test
#'
#' Computes asymptotic critical values and p-values for the modified EM test
#' of Kasahara and Shimotsu for multivariate normal mixtures.
#'
#' @export
#' @title mvnmixCrit
#' @name mvnmixCrit
#' @param y n by d matrix of data
#' @param parlist The parameter estimates as a list containing alpha, mu, and sigma
#'   in the form of (alpha = (alpha_1,...,alpha_m), mu = (mu_1',...,mu_m'),
#'   sigma = (vech(sigma_1)',...,vech(sigma_m)')
#' @param values Vector of test statistic values at which p-values are computed
#' @param nrep Number of simulation replications (default 10000)
#' @param ninits.crit Number of random starts for cone projection optimization (default 25)
#' @return A list with elements `crit` and `pvals`.
mvnmixCrit <- function(y, parlist, values = NULL, nrep = 10000, ninits.crit = 25) {
  y <- mvn_validate_y(y, allow_univariate = TRUE)

  d <- ncol(y)
  alpha <- parlist$alpha
  m <- length(alpha)

  if (d == 1L && m == 1L) {
    crit <- stats::qchisq(c(0.1, 0.05, 0.01), 2, lower.tail = FALSE)
    pvals <- NULL
    if (!is.null(values)) {
      pvals <- stats::pchisq(values, 2, lower.tail = FALSE)
    }
    return(list(crit = crit, pvals = pvals))
  }

  if (d == 1L) {
    return(mvnmixCrit_simulate(y, parlist, values, nrep, ninits.crit = 1L))
  }

  mvnmixCrit_simulate(y, parlist, values, nrep, ninits.crit)
}

mvnmixCrit_simulate <- function(y, parlist, values, nrep, ninits.crit = 25) {
  set.seed(123456)

  n <- nrow(y)
  d <- ncol(y)
  dsig <- d * (d + 1) / 2

  alpha <- parlist$alpha
  mu <- parlist$mu
  sigma <- parlist$sigma
  m <- length(alpha)

  tup3 <- mvn_tuples3(d)
  tup4 <- mvn_tuples4(d)
  d_muv <- nrow(tup3)
  d_mu4 <- nrow(tup4)
  n_lam <- d_muv + d_mu4

  perm12_list <- vector("list", d_muv)
  for (s in seq_len(d_muv)) {
    perm12_list[[s]] <- perm12(tup3[s, 1], tup3[s, 2], tup3[s, 3])
  }

  perm22_list <- vector("list", d_mu4)
  mc4_vec <- numeric(d_mu4)
  for (s in seq_len(d_mu4)) {
    perm22_list[[s]] <- perm22(tup4[s, 1], tup4[s, 2], tup4[s, 3], tup4[s, 4])
    mc4_vec[s] <- perm4_count(tup4[s, 1], tup4[s, 2], tup4[s, 3], tup4[s, 4])
  }

  mu_mat <- matrix(mu, nrow = d, ncol = m)
  Sigma_list <- vector("list", m)
  for (j in seq_len(m)) {
    idx <- ((j - 1L) * dsig + 1L):(j * dsig)
    Sigma_list[[j]] <- sigmavec2mat(sigma[idx], d)
  }

  scores <- mvn_scores(y, alpha, mu_mat, Sigma_list, m, d, tup3, tup4)
  S_eta <- scores$S_eta
  S_lam <- scores$S_lam

  I_eta <- crossprod(S_eta) / n
  I_lam <- crossprod(S_lam) / n
  I_el <- crossprod(S_eta, S_lam) / n

  Iall <- rbind(cbind(I_eta, I_el), cbind(t(I_el), I_lam))
  if (!is.finite(rcond(Iall)) || rcond(Iall) < .Machine$double.eps) {
    Iall <- mvn_regularize_psd(Iall)
    dim_eta <- ncol(S_eta)
    I_eta <- Iall[1:dim_eta, 1:dim_eta, drop = FALSE]
    I_el <- Iall[1:dim_eta, (dim_eta + 1):ncol(Iall), drop = FALSE]
    I_lam <- Iall[(dim_eta + 1):nrow(Iall), (dim_eta + 1):ncol(Iall), drop = FALSE]
  }

  I_eta_inv_Iel <- mvn_stable_solve(I_eta, I_el)
  I_lam_eta <- I_lam - crossprod(I_el, I_eta_inv_Iel)
  I_lam_eta <- mvn_regularize_psd(I_lam_eta)

  eig <- eigen(I_lam_eta, symmetric = TRUE)
  eig$values <- pmax(eig$values, 0)
  sqrtI <- eig$vectors %*% (sqrt(eig$values) * t(eig$vectors))

  u <- matrix(stats::rnorm(nrep * m * n_lam), nrep, m * n_lam) %*% sqrtI
  EM <- matrix(0, nrow = nrep, ncol = m)

  if (d == 1L) {
    for (jj in seq_len(m)) {
      idx <- ((jj - 1L) * n_lam + 1L):(jj * n_lam)
      I_jj <- I_lam_eta[idx, idx, drop = FALSE]
      I_jj_inv <- mvn_stable_solve(I_jj)
      u_jj <- u[, idx, drop = FALSE]
      EM[, jj] <- rowSums((u_jj %*% I_jj_inv) * u_jj)
    }
  } else {
    perm12_flat <- do.call(rbind, perm12_list) - 1L
    perm12_offsets <- c(0L, cumsum(vapply(perm12_list, nrow, integer(1))))

    perm22_flat <- do.call(rbind, perm22_list) - 1L
    perm22_offsets <- c(0L, cumsum(vapply(perm22_list, nrow, integer(1))))

    init_grid <- make_init_grid(d, dsig, ninits.crit)

    EM <- cppConeProjectBatch(
      u = u,
      I_lam_eta = I_lam_eta,
      m = m,
      d = d,
      perm12_flat = perm12_flat,
      perm12_offsets = perm12_offsets,
      perm22_flat = perm22_flat,
      perm22_offsets = perm22_offsets,
      mc4_vec = mc4_vec,
      tup4 = tup4 - 1L,
      init_grid = init_grid,
      d_muv = d_muv,
      d_mu4 = d_mu4
    )
  }

  max_EM <- apply(EM, 1, max)
  max_EM_sort <- sort(max_EM)

  q <- pmin(ceiling(nrep * c(0.90, 0.95, 0.99)), nrep)
  crit <- max_EM_sort[q]

  pvals <- NULL
  if (!is.null(values)) {
    pvals <- colMeans(outer(max_EM_sort, values, ">"))
  }

  list(crit = crit, pvals = pvals)
}

mvn_tuples3 <- function(d) {
  n3 <- choose(d + 2L, 3L)
  out <- matrix(0L, nrow = n3, ncol = 3L)
  idx <- 1L
  for (i in seq_len(d)) {
    for (j in i:d) {
      for (k in j:d) {
        out[idx, ] <- c(i, j, k)
        idx <- idx + 1L
      }
    }
  }
  out
}

mvn_tuples4 <- function(d) {
  n4 <- choose(d + 3L, 4L)
  out <- matrix(0L, nrow = n4, ncol = 4L)
  idx <- 1L
  for (i in seq_len(d)) {
    for (j in i:d) {
      for (k in j:d) {
        for (l in k:d) {
          out[idx, ] <- c(i, j, k, l)
          idx <- idx + 1L
        }
      }
    }
  }
  out
}

perm12 <- function(i, j, k) {
  perms <- rbind(c(i, j, k), c(j, i, k), c(k, i, j))
  perms[, 2:3] <- t(apply(perms[, 2:3, drop = FALSE], 1, sort))
  unique(perms)
}

perm22 <- function(i, j, k, l) {
  idx <- c(i, j, k, l)
  all_perms <- .permutations4(idx)
  out <- matrix(0L, nrow = 0, ncol = 4)
  for (r in seq_len(nrow(all_perms))) {
    p <- all_perms[r, ]
    a <- sort(p[1:2])
    b <- sort(p[3:4])
    if (a[1] > b[1] || (a[1] == b[1] && a[2] > b[2])) {
      tmp <- a
      a <- b
      b <- tmp
    }
    out <- rbind(out, c(a, b))
  }
  unique(out)
}

perm4_count <- function(i, j, k, l) {
  tab <- table(c(i, j, k, l))
  factorial(4) / prod(factorial(tab))
}

.permutations4 <- function(x) {
  n <- length(x)
  if (n == 1L) {
    return(matrix(x, nrow = 1L, ncol = 1L))
  }

  out <- matrix(0L, nrow = 0L, ncol = n)
  for (i in seq_len(n)) {
    rest <- .permutations4(x[-i])
    out <- rbind(out, cbind(x[i], rest))
  }
  unique(out)
}

make_init_grid <- function(d, dsig, ninits) {
  n_par <- d + dsig
  grid <- matrix(0, nrow = ninits, ncol = n_par)
  scales <- c(0.1, 0.5, 1.0, 2.0)

  for (i in seq_len(ninits)) {
    grid[i, ] <- stats::rnorm(n_par, sd = scales[(i %% 4L) + 1L])
  }

  grid
}

mvn_scores <- function(y, alpha, mu_mat, Sigma_list, m, d, tup3, tup4) {
  n <- nrow(y)
  dsig <- d * (d + 1L) / 2L
  d_muv <- nrow(tup3)
  d_mu4 <- nrow(tup4)
  n_lam <- d_muv + d_mu4

  P_list <- vector("list", m)
  z_list <- vector("list", m)
  logf_mat <- matrix(0, n, m)

  for (j in seq_len(m)) {
    Sigma_j <- Sigma_list[[j]]
    chol_obj <- mvn_chol_with_jitter(Sigma_j)
    chol_j <- chol_obj$chol
    P_j <- chol2inv(chol_j)

    P_list[[j]] <- P_j

    resid <- sweep(y, 2, mu_mat[, j], "-")
    z_j <- resid %*% P_j
    z_list[[j]] <- z_j

    log_det <- 2 * sum(log(diag(chol_j)))
    mahal <- rowSums(z_j * resid)
    logf_mat[, j] <- -0.5 * (d * log(2 * pi) + log_det + mahal)
  }

  alpha_safe <- pmax(alpha, 1e-16)
  log_weighted <- sweep(logf_mat, 2, log(alpha_safe), "+")
  row_max <- apply(log_weighted, 1, max)
  weighted_exp <- exp(log_weighted - row_max)
  denom <- rowSums(weighted_exp)
  w_mat <- weighted_exp / denom

  if (m >= 2L) {
    f_over_f0 <- sweep(w_mat, 2, alpha_safe, "/")
    S_alpha <- f_over_f0[, 1:(m - 1L), drop = FALSE] - f_over_f0[, m]
  } else {
    S_alpha <- matrix(0, nrow = n, ncol = 0)
  }

  S_mu <- matrix(0, nrow = n, ncol = m * d)
  for (j in seq_len(m)) {
    S_mu[, ((j - 1L) * d + 1L):(j * d)] <- w_mat[, j] * z_list[[j]]
  }

  S_v <- matrix(0, nrow = n, ncol = m * dsig)
  for (j in seq_len(m)) {
    P_j <- P_list[[j]]
    z_j <- z_list[[j]]
    sv_j <- matrix(0, nrow = n, ncol = dsig)
    idx <- 0L

    for (b in seq_len(d)) {
      for (a in b:d) {
        idx <- idx + 1L
        c_ab <- if (a == b) 1 else 2
        sv_j[, idx] <- w_mat[, j] * 0.5 * c_ab * (z_j[, a] * z_j[, b] - P_j[a, b])
      }
    }

    S_v[, ((j - 1L) * dsig + 1L):(j * dsig)] <- sv_j
  }

  S_eta <- cbind(S_alpha, S_mu, S_v)
  S_lam <- matrix(0, nrow = n, ncol = m * n_lam)

  for (j in seq_len(m)) {
    P_j <- P_list[[j]]
    z_j <- z_list[[j]]

    s_muv <- matrix(0, nrow = n, ncol = d_muv)
    for (s in seq_len(d_muv)) {
      ii <- tup3[s, 1]
      jj <- tup3[s, 2]
      kk <- tup3[s, 3]
      s_muv[, s] <- w_mat[, j] * (
        z_j[, ii] * z_j[, jj] * z_j[, kk] -
          P_j[ii, jj] * z_j[, kk] -
          P_j[ii, kk] * z_j[, jj] -
          P_j[jj, kk] * z_j[, ii]
      ) / 6
    }

    s_mu4 <- matrix(0, nrow = n, ncol = d_mu4)
    for (s in seq_len(d_mu4)) {
      ii <- tup4[s, 1]
      jj <- tup4[s, 2]
      kk <- tup4[s, 3]
      ll <- tup4[s, 4]
      s_mu4[, s] <- w_mat[, j] * (
        z_j[, ii] * z_j[, jj] * z_j[, kk] * z_j[, ll] -
          P_j[ii, jj] * z_j[, kk] * z_j[, ll] -
          P_j[ii, kk] * z_j[, jj] * z_j[, ll] -
          P_j[ii, ll] * z_j[, jj] * z_j[, kk] -
          P_j[jj, kk] * z_j[, ii] * z_j[, ll] -
          P_j[jj, ll] * z_j[, ii] * z_j[, kk] -
          P_j[kk, ll] * z_j[, ii] * z_j[, jj] +
          P_j[ii, jj] * P_j[kk, ll] +
          P_j[ii, kk] * P_j[jj, ll] +
          P_j[ii, ll] * P_j[jj, kk]
      ) / 24
    }

    col_start <- (j - 1L) * n_lam + 1L
    S_lam[, col_start:(col_start + d_muv - 1L)] <- s_muv
    S_lam[, (col_start + d_muv):(col_start + n_lam - 1L)] <- s_mu4
  }

  list(S_eta = S_eta, S_lam = S_lam)
}
