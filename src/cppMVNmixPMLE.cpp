#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

namespace {

constexpr double kDetTol = 1e-12;
constexpr double kAlphaTol = 1e-8;
constexpr double kWeightTol = 1e-12;

inline void vech_to_sym(const arma::vec& src, int offset, int d, arma::mat& out) {
  out.zeros();
  int idx = offset;
  for (int col = 0; col < d; ++col) {
    for (int row = col; row < d; ++row) {
      out(row, col) = src(idx++);
    }
  }
  out = arma::symmatl(out);
}

inline void sym_to_vech(const arma::mat& src, arma::vec& out, int offset) {
  int idx = offset;
  int d = src.n_rows;
  for (int col = 0; col < d; ++col) {
    for (int row = col; row < d; ++row) {
      out(idx++) = src(row, col);
    }
  }
}

} // namespace

// [[Rcpp::export]]
List cppMVNmixPMLE(NumericMatrix bs,
                   NumericMatrix ys,
                   NumericVector mu0s,
                   NumericVector sigma0s,
                   int m,
                   double an,
                   int maxit = 2000,
                   int ninits = 10,
                   double tol = 1e-8,
                   double tau = 0.5,
                   int h = 0,
                   int k = 0) {
  const int n = ys.nrow();
  const int d = ys.ncol();
  const int dsig = d * (d + 1) / 2;

  if (n <= 1 || d <= 0) {
    stop("ys must have at least two rows and one column.");
  }
  if (m <= 0) {
    stop("m must be positive.");
  }
  if (ninits <= 0) {
    stop("ninits must be positive.");
  }
  if (bs.ncol() < ninits) {
    stop("bs must have at least ninits columns.");
  }
  if (bs.nrow() != m + m * d + m * dsig) {
    stop("bs has incompatible row dimension.");
  }
  if (sigma0s.size() != m * dsig) {
    stop("sigma0s has incompatible length.");
  }
  if (k > 0 && (h < 1 || h > m)) {
    stop("h must be in {1, ..., m} when k > 0.");
  }
  if (k == 1 && mu0s.size() < (m + 1)) {
    stop("mu0s must have length at least m + 1 when k == 1.");
  }

  arma::mat b = Rcpp::as<arma::mat>(bs);
  const arma::mat y = Rcpp::as<arma::mat>(ys);
  arma::vec mu0 = Rcpp::as<arma::vec>(mu0s);
  arma::vec sigma0 = Rcpp::as<arma::vec>(sigma0s);

  arma::cube sigma0_cube(d, d, m, arma::fill::zeros);
  arma::cube sigma_cube(d, d, m, arma::fill::zeros);

  arma::vec lb(m, arma::fill::zeros), ub(m, arma::fill::zeros);
  if (k == 1) {
    mu0(0) = R_NegInf;
    mu0(m) = R_PosInf;
    for (int j = 0; j < h; ++j) {
      lb(j) = (mu0(j) + mu0(j + 1)) / 2.0;
      ub(j) = (mu0(j + 1) + mu0(j + 2)) / 2.0;
    }
    for (int j = h; j < m; ++j) {
      lb(j) = (mu0(j - 1) + mu0(j)) / 2.0;
      ub(j) = (mu0(j) + mu0(j + 1)) / 2.0;
    }
  }

  for (int j = 0; j < m; ++j) {
    arma::mat sigma0_j(d, d, arma::fill::zeros);
    vech_to_sym(sigma0, j * dsig, d, sigma0_j);
    sigma0_cube.slice(j) = sigma0_j;
  }

  arma::mat post(m * n, ninits, arma::fill::zeros);
  arma::vec notcg(ninits, arma::fill::zeros);
  arma::vec penloglikset(ninits, arma::fill::value(R_NegInf));
  arma::vec loglikset(ninits, arma::fill::value(R_NegInf));

  arma::vec b_jn(bs.nrow());
  arma::vec alpha(m), mu(m * d), sigma(m * dsig), alp_sig(m);
  arma::vec logdetsigma(m), pen(m);
  arma::mat l_m(n, m, arma::fill::zeros);
  arma::mat r(n, m, arma::fill::zeros);
  arma::vec minr(n, arma::fill::zeros), sum_l_m(n, arma::fill::zeros);
  arma::mat w(n, m, arma::fill::zeros);
  arma::vec mu_j(d), wcol(n);
  arma::mat ydot(n, d), ssr_j(d, d), sigma_j_inv(d, d), s0j(d, d);

  double ll = R_NegInf;
  double penloglik = R_NegInf;

  for (int jn = 0; jn < ninits; ++jn) {
    b_jn = b.col(jn);
    alpha = b_jn.subvec(0, m - 1);
    mu = b_jn.subvec(m, m + m * d - 1);
    sigma = b_jn.subvec(m + m * d, m + m * d + m * dsig - 1);

    for (int j = 0; j < m; ++j) {
      arma::mat sigma_j(d, d, arma::fill::zeros);
      vech_to_sym(sigma, j * dsig, d, sigma_j);
      sigma_cube.slice(j) = sigma_j;
    }

    double oldpenloglik = R_NegInf;
    ll = R_NegInf;
    penloglik = R_NegInf;
    w.zeros();

    bool singular = false;

    for (int iter = 0; iter < maxit; ++iter) {
      for (int j = 0; j < m; ++j) {
        mu_j = mu.subvec(j * d, (j + 1) * d - 1);
        ydot = y.each_row() - mu_j.t();

        const arma::mat& sigma_j = sigma_cube.slice(j);
        bool inv_ok = arma::inv_sympd(sigma_j_inv, sigma_j);
        if (!inv_ok || !sigma_j_inv.is_finite()) {
          singular = true;
          break;
        }

        double logdet_j = 0.0;
        double sign_j = 1.0;
        arma::log_det(logdet_j, sign_j, sigma_j);
        if (sign_j <= 0 || !std::isfinite(logdet_j) || logdet_j < std::log(kDetTol)) {
          singular = true;
          break;
        }

        logdetsigma(j) = logdet_j;

        arma::mat rtilde = 0.5 * (ydot * sigma_j_inv) % ydot;
        r.col(j) = arma::sum(rtilde, 1);

        s0j = sigma0_cube.slice(j) * sigma_j_inv;
        double logdet_s0j = 0.0;
        double sign_s0j = 1.0;
        arma::log_det(logdet_s0j, sign_s0j, s0j);
        if (sign_s0j <= 0 || !std::isfinite(logdet_s0j)) {
          singular = true;
          break;
        }

        pen(j) = arma::trace(s0j) - logdet_s0j - d;
      }

      if (singular) {
        penloglik = R_NegInf;
        ll = R_NegInf;
        break;
      }

      alp_sig = alpha % arma::exp(-0.5 * logdetsigma);
      if (!alp_sig.is_finite() || arma::any(alp_sig <= 0)) {
        singular = true;
        penloglik = R_NegInf;
        ll = R_NegInf;
        break;
      }

      minr = arma::min(r, 1);
      l_m = arma::exp(-(r.each_col() - minr));
      l_m.each_row() %= alp_sig.t();
      sum_l_m = arma::sum(l_m, 1);

      if (!sum_l_m.is_finite() || arma::any(sum_l_m <= kWeightTol)) {
        singular = true;
        penloglik = R_NegInf;
        ll = R_NegInf;
        break;
      }

      w = l_m.each_col() / sum_l_m;

      ll = arma::sum(arma::log(sum_l_m) - minr) - static_cast<double>(n) * d * M_LN_SQRT_2PI;
      penloglik = ll + std::log(2.0) + std::fmin(std::log(tau), std::log(1.0 - tau)) - an * arma::sum(pen);

      double diff = penloglik - oldpenloglik;
      oldpenloglik = penloglik;

      if (diff < tol) {
        break;
      }

      for (int j = 0; j < m; ++j) {
        wcol = w.col(j);
        double w_j = arma::accu(wcol);
        if (!std::isfinite(w_j) || w_j <= kWeightTol) {
          singular = true;
          break;
        }

        alpha(j) = w_j / n;

        mu_j = (y.t() * wcol) / w_j;
        if (k == 1) {
          mu_j(0) = std::fmin(std::fmax(mu_j(0), lb(j)), ub(j));
        }
        mu.subvec(j * d, (j + 1) * d - 1) = mu_j;

        ydot = y.each_row() - mu_j.t();
        ssr_j = (ydot.each_col() % wcol).t() * ydot;

        arma::mat sigma_j = (2.0 * an * sigma0_cube.slice(j) + ssr_j) / (2.0 * an + w_j);
        sigma_cube.slice(j) = sigma_j;
      }

      if (singular) {
        penloglik = R_NegInf;
        ll = R_NegInf;
        break;
      }

      if (k == 1) {
        double alphah = alpha(h - 1) + alpha(h);
        alpha(h - 1) = alphah * tau;
        alpha(h) = alphah * (1.0 - tau);
      } else if (k > 1) {
        double alphah = alpha(h - 1) + alpha(h);
        double tauhat = alpha(h - 1) / alphah;
        if (tauhat <= 0.5) {
          tau = std::fmin((alpha(h - 1) * n + 1.0) / (alpha(h - 1) * n + alpha(h) * n + 1.0), 0.5);
        } else {
          tau = std::fmax(alpha(h - 1) * n / (alpha(h - 1) * n + alpha(h) * n + 1.0), 0.5);
        }
        alpha(h - 1) = alphah * tau;
        alpha(h) = alphah * (1.0 - tau);
      }

      if (!alpha.is_finite() || arma::any(alpha < kAlphaTol)) {
        singular = true;
        penloglik = R_NegInf;
        ll = R_NegInf;
        break;
      }

      for (int j = 0; j < m; ++j) {
        double logdet_j = 0.0;
        double sign_j = 1.0;
        arma::log_det(logdet_j, sign_j, sigma_cube.slice(j));
        if (sign_j <= 0 || !std::isfinite(logdet_j) || logdet_j < std::log(kDetTol)) {
          singular = true;
          penloglik = R_NegInf;
          ll = R_NegInf;
          break;
        }
      }

      if (singular) {
        break;
      }

      for (int j = 0; j < m; ++j) {
        sym_to_vech(sigma_cube.slice(j), sigma, j * dsig);
      }
    }

    if (singular) {
      notcg(jn) = 1;
    }

    penloglikset(jn) = penloglik;
    loglikset(jn) = ll;

    for (int j = 0; j < m; ++j) {
      sym_to_vech(sigma_cube.slice(j), sigma, j * dsig);
    }

    b_jn.subvec(0, m - 1) = alpha;
    b_jn.subvec(m, m + m * d - 1) = mu;
    b_jn.subvec(m + m * d, m + m * d + m * dsig - 1) = sigma;
    b.col(jn) = b_jn;

    post.col(jn) = arma::vectorise(w);
  }

  return List::create(
    Named("penloglikset") = wrap(penloglikset),
    Named("loglikset") = wrap(loglikset),
    Named("notcg") = wrap(notcg),
    Named("post") = wrap(post),
    Named("b") = wrap(b)
  );
}
