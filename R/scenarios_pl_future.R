#' scenarios function (fully optimized with Rcpp)
#' This function computes the mean and covariances to draw from the conditional forecast
#' The actual draw is done in the simscen function
#'
#' @param h forecast horizon
#' @param path conditional path of observables
#' @param obs position of observable(s)
#' @param free_shocks position of non-driving shocks (NA if all driving)
#' @param n_sample Number of draws to sample (<= n_draws)
#' @param n_cores Number of parallel workers (<= n_draws)
#' @param data_ Optional matrix of data (default taken from matrices$Z). Note: the last observation
#'        in data_ is used as the starting point; it should not overlap with the scenario forecasting period.
#' @param g Optional matrix of non-driving shocks
#' @param Sigma_g Optional covariance matrix of non-driving shocks
#' @param posterior Optional posterior object (default taken from calling environment)
#' @param matrices Optional matrices object from gen_mats() (default taken from calling environment)
#'
#' @return list of mu_eps, Sigma_eps, mu_y, Sigma_y, big_b, big_M, draws_used
#' @examples
#' \donttest{
#' library(APRScenario)
#' data(NKdata)
#'
#' # Minimal example with a toy specification
#' spec <- bsvarSIGNs::specify_bsvarSIGN$new(as.matrix(NKdata[,2:4]), p = 1)
#' posterior <- bsvars::estimate(spec, S = 10)  # Use small S for fast test
#' matrices <- gen_mats(posterior = posterior, specification = spec)
#' # and having posterior object
#' scenario_result <- scenarios(h = 2,
#'                              path = c(1.0, 1.1),
#'                              obs = 1,
#'                              free_shocks = NA,
#'                              posterior = posterior,
#'                              matrices = matrices)
#' }
#' @export
#' @useDynLib APRScenario, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppParallel

scenarios_pl_future <- function(h = 3,
                                path = NULL,
                                obs = NULL,
                                free_shocks = NULL,
                                n_sample = NULL,
                                n_cores = NULL,
                                data_ = NULL,
                                g = NULL, Sigma_g = NULL,
                                posterior = NULL,
                                matrices = NULL) {
  
  # Get matrices from calling environment if not provided
  if (is.null(matrices)) {
    if (exists("matrices", envir = parent.frame())) {
      matrices <- get("matrices", envir = parent.frame())
    } else {
      stop("Please provide matrices object from gen_mats() or ensure it exists in calling environment")
    }
  }
  
  # Get posterior from calling environment if not provided
  if (is.null(posterior)) {
    if (exists("posterior", envir = parent.frame())) {
      posterior <- get("posterior", envir = parent.frame())
    } else {
      stop("Please provide posterior object or ensure it exists in calling environment")
    }
  }
  
  # Get data from matrices if not provided
  if (is.null(data_)) {
    data_ <- matrices$Z
  }
  
  stopifnot(length(path) == length(obs) * h)
  if (is.null(dim(path))) {
    if (length(path) != h) {
      stop("Length of path (", length(path), ") is not equal to h (", h, ")")
    }
  } else {
    if (dim(path)[1] != length(obs)) {
      stop("path must be n_constrained_vars x h, got ", dim(path)[1], " x ", dim(path)[2])
    }
  }
  
  n_var  <- dim(posterior$posterior$B)[1]
  n_p    <- (dim(posterior$posterior$A)[2] - 1) / n_var
  n_draws <- dim(posterior$posterior$B)[3]
  if (is.null(n_sample)) n_sample <- n_draws
  
  # Build big_b, big_M as before (single-threaded C++ inside)
  tmp   <- big_b_and_M_future(h = h, n_draws = n_draws, n_var = n_var,
                          n_p = n_p, data_ = data_, matrices = matrices, n_cores = n_cores)
  big_b <- tmp[[1]]
  big_M <- tmp[[2]]
  
  draws_to_use <- if (n_sample < n_draws) sample(seq_len(n_draws), n_sample) else seq_len(n_draws)
  big_b        <- big_b[, , draws_to_use, drop = FALSE]
  big_M        <- big_M[, , draws_to_use, drop = FALSE]
  n_draws      <- n_sample
  
  shock_idx <- if (any(is.na(free_shocks))) NA_integer_ else as.integer(free_shocks)
  if (!is.null(g)) g <- as.numeric(t(g)) # flatten g
  
  # Set up future plan: multisession on Windows, workers = n_cores or availableCores()
  if (is.null(n_cores)) {
    n_cores <- future::availableCores()
  }
  future::plan(future::multisession, workers = n_cores)
  
  # One-draw worker: passes one slice to C++ core
  one_draw_fun <- function(i) {
    b_i <- big_b[, , i]
    M_i <- big_M[, , i]
    full_scenario_core_one(b = b_i,
                           M = M_i,
                           obs = as.integer(obs),
                           path = as.numeric(path),
                           shocks = shock_idx,
                           h = h,
                           n_var = n_var,
                           g_ = g,
                           Sigma_g_ = Sigma_g)
  }
  
  # Parallel over draws
  out_list <- future.apply::future_lapply(seq_len(n_draws), one_draw_fun)
  
  nM <- dim(big_M)[1]
  
  # Extract and bind arrays
  mu_eps_list   <- lapply(out_list, `[[`, "mu_eps")
  Sigma_eps_list<- lapply(out_list, `[[`, "Sigma_eps")
  mu_y_list     <- lapply(out_list, `[[`, "mu_y")
  Sigma_y_list  <- lapply(out_list, `[[`, "Sigma_y")
  
  list(
    mu_eps   = abind::abind(lapply(mu_eps_list,   function(x) matrix(x, ncol = 1)), along = 3),
    Sigma_eps= abind::abind(lapply(Sigma_eps_list,function(x) matrix(x, nrow = nM, ncol = nM)), along = 3),
    mu_y    = abind::abind(lapply(mu_y_list,     function(x) matrix(x, ncol = 1)), along = 3),
    Sigma_y = abind::abind(lapply(Sigma_y_list,  function(x) matrix(x, nrow = nM, ncol = nM)), along = 3),
    big_b   = big_b,
    big_M   = big_M,
    draws_used = draws_to_use
  )
}