#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// ============================================================
// vech -> symmetric matrix  (lower-tri, column-major, same as R)
// ============================================================
static mat vech_to_symmat(const double* v, int d) {
  mat S(d, d, fill::zeros);
  int idx = 0;
  for (int col = 0; col < d; col++) {
    for (int row = col; row < d; row++) {
      S(row, col) = v[idx++];
    }
  }
  S = S + S.t();
  S.diag() *= 0.5;
  return S;
}

// ============================================================
// Cone 1 map: (lam_mu, lam_v) -> t_vec  of length d_muv + d_mu4
//   t_muv[s] = sum over perms12[s]: lam_mu[p1] * lam_v_mat(p2,p3)
//   t_mu4[s] = sum over perms22[s]: lam_v_mat(p1,p2) * lam_v_mat(p3,p4)
// ============================================================
static vec cone1_map(const vec& par, int d, int dsig, int d_muv, int d_mu4,
                     const imat& perm12_flat, const ivec& perm12_off,
                     const imat& perm22_flat, const ivec& perm22_off) {
  vec lam_mu = par.head(d);
  mat lam_v = vech_to_symmat(par.memptr() + d, d);

  vec out(d_muv + d_mu4);

  // mu-v part
  for (int s = 0; s < d_muv; s++) {
    double val = 0.0;
    int r0 = perm12_off(s), r1 = perm12_off(s + 1);
    for (int r = r0; r < r1; r++) {
      val += lam_mu(perm12_flat(r, 0)) * lam_v(perm12_flat(r, 1), perm12_flat(r, 2));
    }
    out(s) = val;
  }

  // mu4 part
  for (int s = 0; s < d_mu4; s++) {
    double val = 0.0;
    int r0 = perm22_off(s), r1 = perm22_off(s + 1);
    for (int r = r0; r < r1; r++) {
      val += lam_v(perm22_flat(r, 0), perm22_flat(r, 1)) *
             lam_v(perm22_flat(r, 2), perm22_flat(r, 3));
    }
    out(d_muv + s) = val;
  }

  return out;
}

// ============================================================
// Cone 2 map: (lam_mu, lam_v) -> t_vec
//   t_muv same as cone1
//   t_mu4[s] = -mc4[s] * prod(lam_mu[tup4[s,]])
// ============================================================
static vec cone2_map(const vec& par, int d, int dsig, int d_muv, int d_mu4,
                     const imat& perm12_flat, const ivec& perm12_off,
                     const vec& mc4_vec, const imat& tup4) {
  vec lam_mu = par.head(d);
  mat lam_v = vech_to_symmat(par.memptr() + d, d);

  vec out(d_muv + d_mu4);

  // mu-v part (identical to cone1)
  for (int s = 0; s < d_muv; s++) {
    double val = 0.0;
    int r0 = perm12_off(s), r1 = perm12_off(s + 1);
    for (int r = r0; r < r1; r++) {
      val += lam_mu(perm12_flat(r, 0)) * lam_v(perm12_flat(r, 1), perm12_flat(r, 2));
    }
    out(s) = val;
  }

  // mu4 part
  for (int s = 0; s < d_mu4; s++) {
    out(d_muv + s) = -mc4_vec(s) *
      lam_mu(tup4(s, 0)) * lam_mu(tup4(s, 1)) *
      lam_mu(tup4(s, 2)) * lam_mu(tup4(s, 3));
  }

  return out;
}

// ============================================================
// Nelder-Mead (minimization, no constraints)
// Standard coefficients: alpha=1, gamma=2, rho=0.5, sigma=0.5
// ============================================================

// Objective function type
typedef double (*ObjFn)(const vec& par, const void* data);

// Cone objective data bundle
struct ConeObjData {
  const vec* Z;
  const mat* W;
  int d, dsig, d_muv, d_mu4;
  const imat* perm12_flat;
  const ivec* perm12_off;
  const imat* perm22_flat;
  const ivec* perm22_off;
  const vec* mc4_vec;
  const imat* tup4;
  int cone_id; // 1 or 2
};

static double cone_obj(const vec& par, const void* data) {
  const ConeObjData* D = static_cast<const ConeObjData*>(data);
  vec t_vec;
  if (D->cone_id == 1) {
    t_vec = cone1_map(par, D->d, D->dsig, D->d_muv, D->d_mu4,
                      *D->perm12_flat, *D->perm12_off,
                      *D->perm22_flat, *D->perm22_off);
  } else {
    t_vec = cone2_map(par, D->d, D->dsig, D->d_muv, D->d_mu4,
                      *D->perm12_flat, *D->perm12_off,
                      *D->mc4_vec, *D->tup4);
  }
  vec diff = t_vec - *D->Z;
  return dot(diff, *D->W * diff);
}

static double nelder_mead(ObjFn fn, const vec& x0, const void* data,
                          int maxiter = 200, double tol = 1e-10) {
  const double alpha_nm = 1.0, gamma_nm = 2.0, rho_nm = 0.5, sigma_nm = 0.5;
  int n = (int)x0.n_elem;
  int np1 = n + 1;

  // Initialize simplex with scale-appropriate steps
  mat simplex(n, np1);
  vec fvals(np1);
  simplex.col(0) = x0;
  fvals(0) = fn(x0, data);
  for (int i = 0; i < n; i++) {
    vec xi = x0;
    double step = (std::abs(x0(i)) > 0.1) ? 0.5 * x0(i) : 0.5;
    xi(i) += step;
    simplex.col(i + 1) = xi;
    fvals(i + 1) = fn(xi, data);
  }

  for (int iter = 0; iter < maxiter; iter++) {
    // Sort
    uvec order = sort_index(fvals);
    mat s_sorted(n, np1);
    vec f_sorted(np1);
    for (int i = 0; i < np1; i++) {
      s_sorted.col(i) = simplex.col(order(i));
      f_sorted(i) = fvals(order(i));
    }
    simplex = s_sorted;
    fvals = f_sorted;

    // Check convergence
    double frange = fvals(np1 - 1) - fvals(0);
    if (frange < tol && frange >= 0) break;

    // Centroid (exclude worst)
    vec centroid = mean(simplex.cols(0, n - 1), 1);
    vec worst = simplex.col(n);

    // Reflection
    vec xr = centroid + alpha_nm * (centroid - worst);
    double fr = fn(xr, data);

    if (fr < fvals(0)) {
      // Expansion
      vec xe = centroid + gamma_nm * (xr - centroid);
      double fe = fn(xe, data);
      if (fe < fr) {
        simplex.col(n) = xe;
        fvals(n) = fe;
      } else {
        simplex.col(n) = xr;
        fvals(n) = fr;
      }
    } else if (fr < fvals(n - 1)) {
      simplex.col(n) = xr;
      fvals(n) = fr;
    } else {
      // Contraction
      vec xc;
      double fc;
      if (fr < fvals(n)) {
        // Outside contraction
        xc = centroid + rho_nm * (xr - centroid);
        fc = fn(xc, data);
      } else {
        // Inside contraction
        xc = centroid + rho_nm * (worst - centroid);
        fc = fn(xc, data);
      }
      if (fc < std::min(fr, fvals(n))) {
        simplex.col(n) = xc;
        fvals(n) = fc;
      } else {
        // Shrink
        vec best = simplex.col(0);
        for (int i = 1; i < np1; i++) {
          simplex.col(i) = best + sigma_nm * (simplex.col(i) - best);
          fvals(i) = fn(simplex.col(i), data);
        }
      }
    }
  }

  // Return best value and optionally the best point
  uvec final_order = sort_index(fvals);
  return fvals(final_order(0));
}

// Nelder-Mead returning both value and best point
static double nelder_mead_xbest(ObjFn fn, const vec& x0, const void* data,
                                vec& xbest, int maxiter = 200, double tol = 1e-10) {
  const double alpha_nm = 1.0, gamma_nm = 2.0, rho_nm = 0.5, sigma_nm = 0.5;
  int n = (int)x0.n_elem;
  int np1 = n + 1;

  mat simplex(n, np1);
  vec fvals(np1);
  simplex.col(0) = x0;
  fvals(0) = fn(x0, data);
  for (int i = 0; i < n; i++) {
    vec xi = x0;
    double step = (std::abs(x0(i)) > 0.1) ? 0.5 * x0(i) : 0.5;
    xi(i) += step;
    simplex.col(i + 1) = xi;
    fvals(i + 1) = fn(xi, data);
  }

  for (int iter = 0; iter < maxiter; iter++) {
    uvec order = sort_index(fvals);
    mat s_sorted(n, np1);
    vec f_sorted(np1);
    for (int i = 0; i < np1; i++) {
      s_sorted.col(i) = simplex.col(order(i));
      f_sorted(i) = fvals(order(i));
    }
    simplex = s_sorted;
    fvals = f_sorted;

    double frange = fvals(np1 - 1) - fvals(0);
    if (frange < tol && frange >= 0) break;

    vec centroid = mean(simplex.cols(0, n - 1), 1);
    vec worst = simplex.col(n);
    vec xr = centroid + alpha_nm * (centroid - worst);
    double fr = fn(xr, data);

    if (fr < fvals(0)) {
      vec xe = centroid + gamma_nm * (xr - centroid);
      double fe = fn(xe, data);
      if (fe < fr) { simplex.col(n) = xe; fvals(n) = fe; }
      else { simplex.col(n) = xr; fvals(n) = fr; }
    } else if (fr < fvals(n - 1)) {
      simplex.col(n) = xr; fvals(n) = fr;
    } else {
      vec xc; double fc;
      if (fr < fvals(n)) {
        xc = centroid + rho_nm * (xr - centroid); fc = fn(xc, data);
      } else {
        xc = centroid + rho_nm * (worst - centroid); fc = fn(xc, data);
      }
      if (fc < std::min(fr, fvals(n))) {
        simplex.col(n) = xc; fvals(n) = fc;
      } else {
        vec best = simplex.col(0);
        for (int i = 1; i < np1; i++) {
          simplex.col(i) = best + sigma_nm * (simplex.col(i) - best);
          fvals(i) = fn(simplex.col(i), data);
        }
      }
    }
  }

  uvec final_order = sort_index(fvals);
  xbest = simplex.col(final_order(0));
  return fvals(final_order(0));
}

// ============================================================
// Single cone projection for one draw
// ============================================================
static double cone_project_single(const vec& Z, const mat& W, double LR_full,
                                  int d, int dsig, int d_muv, int d_mu4,
                                  const imat& perm12_flat, const ivec& perm12_off,
                                  const imat& perm22_flat, const ivec& perm22_off,
                                  const vec& mc4_vec, const imat& tup4,
                                  const mat& init_grid) {
  int n_par = d + dsig;
  int ninits = (int)init_grid.n_rows;

  ConeObjData data;
  data.Z = &Z;
  data.W = &W;
  data.d = d;
  data.dsig = dsig;
  data.d_muv = d_muv;
  data.d_mu4 = d_mu4;
  data.perm12_flat = &perm12_flat;
  data.perm12_off = &perm12_off;
  data.perm22_flat = &perm22_flat;
  data.perm22_off = &perm22_off;
  data.mc4_vec = &mc4_vec;
  data.tup4 = &tup4;

  vec zero_par(n_par, fill::zeros);

  // Cone 1: multi-start NM + restart from best
  data.cone_id = 1;
  double best1 = cone_obj(zero_par, &data);
  vec best1_x = zero_par;
  for (int i = 0; i < ninits; i++) {
    vec par0 = init_grid.row(i).t();
    vec xb(n_par);
    double val = nelder_mead_xbest(cone_obj, par0, &data, xb);
    if (val < best1) { best1 = val; best1_x = xb; }
  }
  // Restart from best point found
  best1 = nelder_mead(cone_obj, best1_x, &data, 200, 1e-10);

  // Cone 2: multi-start NM + restart from best
  data.cone_id = 2;
  double best2 = cone_obj(zero_par, &data);
  vec best2_x = zero_par;
  for (int i = 0; i < ninits; i++) {
    vec par0 = init_grid.row(i).t();
    vec xb(n_par);
    double val = nelder_mead_xbest(cone_obj, par0, &data, xb);
    if (val < best2) { best2 = val; best2_x = xb; }
  }
  // Restart from best point found
  best2 = nelder_mead(cone_obj, best2_x, &data, 200, 1e-10);

  double LR1 = std::max(LR_full - best1, 0.0);
  double LR2 = std::max(LR_full - best2, 0.0);
  return std::max(LR1, LR2);
}

// ============================================================
// Exported batch function
// ============================================================
// [[Rcpp::export]]
arma::mat cppConeProjectBatch(arma::mat u,
                              arma::mat I_lam_eta,
                              int m, int d,
                              arma::imat perm12_flat,
                              arma::ivec perm12_offsets,
                              arma::imat perm22_flat,
                              arma::ivec perm22_offsets,
                              arma::vec mc4_vec,
                              arma::imat tup4,
                              arma::mat init_grid,
                              int d_muv, int d_mu4) {
  int nrep = (int)u.n_rows;
  int dsig = d * (d + 1) / 2;
  int n_lam = d_muv + d_mu4;

  mat EM(nrep, m, fill::zeros);

  for (int jj = 0; jj < m; jj++) {
    int idx0 = jj * n_lam;
    int idx1 = idx0 + n_lam - 1;

    mat I_jj = I_lam_eta.submat(idx0, idx0, idx1, idx1);
    mat I_jj_inv;
    bool ok = inv_sympd(I_jj_inv, I_jj);
    if (!ok || !I_jj_inv.is_finite()) {
      I_jj_inv = pinv(I_jj);
    }
    mat u_jj = u.cols(idx0, idx1);
    mat Z_jj = u_jj * I_jj_inv;

    // W = I_jj (the information submatrix)
    mat W = I_jj;

    // Cholesky for LR_full computation (not strictly needed; use W directly)
    for (int rr = 0; rr < nrep; rr++) {
      vec Z_r = Z_jj.row(rr).t();
      double LR_full = dot(u_jj.row(rr), Z_jj.row(rr));

      EM(rr, jj) = cone_project_single(Z_r, W, LR_full,
                                        d, dsig, d_muv, d_mu4,
                                        perm12_flat, perm12_offsets,
                                        perm22_flat, perm22_offsets,
                                        mc4_vec, tup4, init_grid);
    }
  }

  return EM;
}
