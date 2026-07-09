#' big_b_and_M
#' This function returns the extended b and M matrices as in APR
#'
#' @param h forecast horison
#' @param n_draws Number of draws
#' @param n_var Number of variables
#' @param n_p Number of lags
#' @param n_core Number of parallel workers
#' @param data_ (matrix optional) The data, stacking Y over X (data and laggs)
#'        -- columns are observations (default taken from matrices$Z)
#'        NB: this is not necessarily the same as the data used to estimate the model
#'        If run counterfactuals in previoius historical period (ie not forecast) must pass the data up to previous period relative to counterfactual
#' @param matrices Optional matrices object from gen_mats() (default taken from calling environment)
#' @returns the big_b and big_M matrices of mean and IRF
#' @examples
#' \dontrun{
#' # Example usage for creating extended matrices
#' result <- big_b_and_M(h = 4, n_draws = 1000, n_var = 3, n_p = 2,
#'                       matrices = matrices)
#' big_b <- result[[1]]
#' big_M <- result[[2]]
#' }
#' @export
#' @import dplyr
#' Parallel big_b_and_M using future_lapply over draws
#' @export
big_b_and_M_future <- function(h, n_draws, n_var, n_p,
                               data_ = NULL, matrices = NULL,
                               n_cores = NULL) {
  
  # Get matrices from calling environment if not provided
  if (is.null(matrices)) {
    if (exists("matrices", envir = parent.frame())) {
      matrices <- get("matrices", envir = parent.frame())
    } else {
      stop("Please provide matrices object from gen_mats() or ensure it exists in calling environment")
    }
  }
  
  # Get data from matrices if not provided
  if (is.null(data_)) {
    data_ <- matrices$Z
  }
  
  # Draws to use (here: 1..n_draws, adapt if you want subsampling)
  draws_to_use <- seq_len(n_draws)
  
  # Set up future plan
  if (is.null(n_cores)) {
    n_cores <- future::availableCores()
  }
  future::plan(future::multisession, workers = n_cores)
  
  # Parallel over draws: each call returns b_h (length n_var) and M_h (list length h)
  per_draw <- future.apply::future_lapply(
    draws_to_use,
    function(d) {
      mat_forc_one(
        h       = h,
        n_var   = n_var,
        n_p     = n_p,
        data_   = data_,
        matrices = matrices,
        d       = d
      )
    }
  )
  
  ## Assemble big_b: 1 x (n_var * h) x n_draws
  # First build b_h_array: 1 x n_var x n_draws at horizon h
  b_h_array <- abind::abind(
    lapply(per_draw, function(x) matrix(x$b_h, nrow = 1)),
    along = 3
  )
  big_b <- array(0, dim = c(1, n_var * h, n_draws))
  
  # Fill big_b block by block over horizons, following your original logic
  for (cnt in seq_len(h)) {
    # For horizon cnt, you want the mean at that horizon; here we use b_h_array (final horizon)
    big_b[1, (1 + n_var * (cnt - 1)):(cnt * n_var), ] <- b_h_array
  }
  
  ## Assemble big_M: (n_var * h) x (n_var * h) x n_draws
  big_M <- array(0, dim = c(n_var * h, n_var * h, n_draws))
  
  # First aggregate M_h over draws for each horizon cnt: n_var x n_var x n_draws
  M_h_global <- vector("list", h)
  for (cnt in seq_len(h)) {
    M_h_global[[cnt]] <- abind::abind(
      lapply(per_draw, function(x) x$M_h[[cnt]]),
      along = 3
    )
  }
  
  # Then fill big_M according to your original block structure
  for (cnt in seq_len(h)) {
    zz <- 1
    for (cnt2 in cnt:h) {
      big_M[(1 + n_var * (cnt - 1)):(cnt * n_var),
            (1 + n_var * (cnt2 - 1)):(cnt2 * n_var), ] <- M_h_global[[zz]]
      zz <- zz + 1
    }
  }
  
  list(big_b = big_b,
       big_M = big_M)
}