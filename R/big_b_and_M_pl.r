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
#' Parallel big_b_and_M using lapply over draws
#' @export
big_b_and_M_pl <- function(h, n_draws, n_var, n_p,
                               data_ = NULL, matrices = NULL,
                               n_cores = NULL) {
  
  if (is.null(matrices)) {
    if (exists("matrices", envir = parent.frame())) {
      matrices <- get("matrices", envir = parent.frame())
    } else {
      stop("Please provide matrices object from gen_mats() or ensure it exists in calling environment")
    }
  }
  
  if (is.null(data_)) {
    data_ <- matrices$Z
  }
  
  if (is.null(n_cores)) {
    n_cores <- max(1L, parallel::detectCores() - 1L)
  }
  
  n_cores <- min(as.integer(n_cores), n_draws)
  draws_to_use <- seq_len(n_draws)
  
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(
    cl,
    varlist = c("h", "n_var", "n_p", "data_", "matrices", "mat_forc_one"),
    envir = environment()
  )
  
  per_draw <- parallel::parLapply(cl, draws_to_use, function(d) {
    out <- mat_forc_one(
      h = h,
      n_var = n_var,
      n_p = n_p,
      data_ = data_,
      matrices = matrices,
      d = d
    )
    
    list(
      b_h = out[[1L]],
      M_h = out[[2L]]
    )
  })
  
  big_b <- abind::abind(
    lapply(per_draw, function(x) {
      array(x[[1L]], dim = c(1, n_var * h, 1))
    }),
    along = 3
  )
  
  big_M <- array(0, dim = c(n_var * h, n_var * h, n_draws))
  
  for (d in seq_len(n_draws)) {
    M_h_draw <- per_draw[[d]][[2L]]    
    for (cnt in seq_len(h)) {
      zz <- 1
      
      for (cnt2 in cnt:h) {
        big_M[
          (1 + n_var * (cnt - 1)):(cnt * n_var),
          (1 + n_var * (cnt2 - 1)):(cnt2 * n_var),
          d
        ] <- M_h_draw[[zz]]
        
        zz <- zz + 1
      }
    }
  }
  
  list(big_b = big_b, big_M = big_M)
}