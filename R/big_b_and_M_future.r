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
    n_cores <- future::availableCores()
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1 || is.na(n_cores) || n_cores < 1) {
    stop("n_cores must be a single integer >= 1")
  }
  n_cores <- as.integer(n_cores)
  
  draws_to_use <- seq_len(n_draws)
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  if (n_cores == 1L) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = n_cores)
  }
  
  per_draw <- future.apply::future_lapply(
    draws_to_use,
    function(d) {
      out <- mat_forc_one(
        h = h,
        n_var = n_var,
        n_p = n_p,
        data_ = data_,
        matrices = matrices,
        d = d
      )
      
      if (is.null(out)) {
        stop(sprintf("mat_forc_one returned NULL for draw %d", d))
      }
      
      if (is.null(out$b_h) && !is.null(out[[1]])) {
        out$b_h <- out[[1]]
      }
      if (is.null(out$M_h) && !is.null(out[[2]])) {
        out$M_h <- out[[2]]
      }
      
      if (is.null(out$b_h)) {
        stop(sprintf("Missing b_h for draw %d", d))
      }
      if (is.null(out$M_h)) {
        stop(sprintf("Missing M_h for draw %d", d))
      }
      if (length(out$b_h) != n_var * h) {
        stop(sprintf("b_h has length %d for draw %d, expected %d", length(out$b_h), d, n_var * h))
      }
      if (!is.list(out$M_h) || length(out$M_h) != h) {
        stop(sprintf("M_h must be a list of length h=%d for draw %d", h, d))
      }
      
      out
    },
    future.seed = TRUE
  )
  
  big_b <- abind::abind(
    lapply(per_draw, function(x) array(x$b_h, dim = c(1, n_var * h, 1))),
    along = 3
  )
  
  big_M <- array(0, dim = c(n_var * h, n_var * h, n_draws))
  
  for (d in seq_len(n_draws)) {
    M_h_draw <- per_draw[[d]]$M_h
    for (cnt in seq_len(h)) {
      zz <- 1
      for (cnt2 in cnt:h) {
        block <- M_h_draw[[zz]]
        if (!is.matrix(block) || !all(dim(block) == c(n_var, n_var))) {
          stop(sprintf(
            "M_h[[%d]] for draw %d must be a %d x %d matrix",
            zz, d, n_var, n_var
          ))
        }
        big_M[(1 + n_var * (cnt - 1)):(cnt * n_var),
              (1 + n_var * (cnt2 - 1)):(cnt2 * n_var),
              d] <- block
        zz <- zz + 1
      }
    }
  }
  
  list(big_b = big_b, big_M = big_M)
}