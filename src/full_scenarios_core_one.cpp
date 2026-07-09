// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List full_scenario_core_one(const arma::mat& b,
                            const arma::mat& M,
                            const IntegerVector& obs,
                            const NumericVector& path,
                            const IntegerVector& shocks,
                            int h,
                            int n_var,
                            Nullable<arma::vec> g_ = R_NilValue,
                            Nullable<arma::mat> Sigma_g_ = R_NilValue) {

  int k_0 = obs.size() * h;
  bool has_shocks = std::any_of(shocks.begin(), shocks.end(), [](int val) {
    return val != NA_INTEGER;
  });
  int k_s = has_shocks ? shocks.size() * h : 0;

  vec g = g_.isNotNull() ? as<vec>(g_) : vec(k_s, fill::zeros);

  mat M_t = trans(M);
  mat b_t = trans(b);

  mat C_h(k_0, n_var * h, fill::zeros);
  for (int j = 0; j < obs.size(); ++j) {
    int var_idx = obs[j] - 1;
    for (int t = 0; t < h; ++t) {
      C_h(j * h + t, var_idx + t * n_var) = 1.0;
    }
  }

  mat f(k_0 + k_s, 1, fill::zeros);
  if (path.size() != obs.size() * h) {
    stop("path must be of length obs × h");
  }
  for (int j = 0; j < obs.size(); ++j) {
    for (int t = 0; t < h; ++t) {
      int i_f = j * h + t;
      int i_path = t * obs.size() + j; // column-major reshape
      f(i_f, 0) = path[i_path];
    }
  }

  mat C_hat = C_h;
  if (has_shocks) {
    mat Xi(k_s, n_var * h, fill::zeros);
    for (int j = 0; j < shocks.size(); ++j) {
      if (shocks[j] == NA_INTEGER) continue;
      int var_idx = shocks[j] - 1;
      for (int t = 0; t < h; ++t) {
        Xi(j * h + t, var_idx + t * n_var) = 1.0;
      }
    }
    mat C_l = Xi * inv(M_t);
    C_hat = join_vert(C_h, C_l);
    f.rows(k_0, k_0 + k_s - 1) = C_l * b_t + g;
  }

  mat D = C_hat * M_t;
  mat D_ast = pinv(D); // Moore–Penrose pseudo-inverse

  mat Omega_f_hat;
  if (Sigma_g_.isNotNull()) {
    mat Sigma_g = as<mat>(Sigma_g_);
    mat Omega_f = Sigma_g;
    mat Z0(k_0, k_s, fill::zeros);
    mat Z1(k_s, k_0, fill::zeros);
    mat I_ks = eye(k_s, k_s);
    Omega_f_hat = join_vert(
      join_horiz(Omega_f, Z0),
      join_horiz(Z1, I_ks)
    );
  } else {
    Omega_f_hat = D * D.t();
  }

  mat mu_eps = D_ast * (f - C_hat * b_t);

  mat I = eye(D_ast.n_rows, D_ast.n_rows);
  mat Sigma_eps = D_ast * Omega_f_hat * D_ast.t()
                  + (I - D_ast * D) * (I - D_ast * D).t();

  vec eps_draw = mu_eps;
  if (has_shocks && Sigma_g_.isNotNull()) {
    mat L_g = chol(as<mat>(Sigma_g_), "lower");
    vec z = randn(k_s);
    vec delta = g + L_g * z;
    for (int j = 0; j < shocks.size(); ++j) {
      if (shocks[j] == NA_INTEGER) continue;
      int var_idx = shocks[j] - 1;
      for (int t = 0; t < h; ++t) {
        int i_eps   = var_idx + t * n_var;
        int i_delta = j * h + t;
        eps_draw(i_eps) = delta(i_delta);
      }
    }
  }

  mat mu_y = b_t + M_t * mu_eps;
  mat Sigma_y = M_t * M + (M_t * D_ast) * (Omega_f_hat - D * D.t()) * (D_ast.t() * M);

  // Optional diagnostic: ensure conditional forecast matches imposed path
  vec mu_y_proj = C_h * mu_y;
  for (arma::uword i = 0; i < mu_y_proj.n_elem; ++i) {
    if (std::abs(mu_y_proj(i) - f(i, 0)) > 1e-6) {
      Rcpp::warning("Mismatch at constraint %d: predicted=%.8f, imposed=%.8f",
                    i + 1, mu_y_proj(i), f(i, 0));
    }
  }

  return List::create(
    Named("mu_eps")   = mu_eps,
    Named("Sigma_eps")= Sigma_eps,
    Named("mu_y")     = mu_y,
    Named("Sigma_y")  = Sigma_y
  );
}