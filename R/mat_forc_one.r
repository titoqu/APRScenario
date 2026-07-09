#' mat_forc_future function
#'##############################################################################
#'  NB: HERE WE USE Antolin-Diaz et al notation                                #
#'  B is reduced form;                                                         #
#'  A is structural;                                                           #
#'  d is intercepts                                                            #
#'  M is reduced so that E(u*u')=Sigma=(A_0*A_0')^(-1) and M_0=A_0^(-1)*Q      #
#'  Note that the code returns conflicting notation:                           #
#'  B=>A_0^(-1)*Q and                                                          #
#'  A=>B                                                                       #
#'##############################################################################
#' @param h Integer; forecast horizon.                                         
#' @param n_draws Integer; number of posterior draws to use.                   
#' @param n_var Integer; number of endogenous variables.                       
#' @param n_p Integer; number of lags in the VAR.                              
#' @param n_cores Integer; number of parallel workers for `future_lapply()`.   
#' @param data_ Optional matrix of data, stacking Y over X 
#' @param matrices Optional matrices object from `gen_mats()` (default taken   
#'        from the calling environment).                                       
#' @returns the big_b and big_M matrices of mean and IRF for a single draw d
#' @examples
#' \donttest{
#' library(APRScenario)
#' data(NKdata)
#'
#' # Minimal example with a toy specification
#' spec <- bsvarSIGNs::specify_bsvarSIGN$new(as.matrix(NKdata[, 2:4]), p = 1)
#' est  <- bsvars::estimate(spec, S = 10)  # Use small S for fast test
#' matrices <- gen_mats(posterior = est, specification = spec)
#'
#' n_var   <- dim(est$posterior$B)[1]
#' n_p     <- (dim(est$posterior$A)[2] - 1) / n_var
#' n_draws <- dim(est$posterior$B)[3]
#'
#' # Future-based construction of big_b and big_M over draws
#' tmp <- big_b_and_M_future(h = 4, n_draws = n_draws,
#'                           n_var = n_var, n_p = n_p,
#'                           matrices = matrices, n_cores = 2)
#' big_b <- tmp$big_b
#' big_M <- tmp$big_M
#' }
#' @export
#' @import dplyr

mat_forc_one <- function(h, n_var, n_p, data_, matrices, d) {
  # Build a "single-draw" matrices object
  matrices_d <- matrices
  
  # Subset M and intercept to draw d
  matrices_d$M         <- matrices$M[, , d, drop = FALSE]
  matrices_d$intercept <- matrices$intercept[, d, drop = FALSE]
  
  # Subset each B_list lag to draw d
  matrices_d$B_list <- lapply(matrices$B_list, function(Bj) {
    Bj[, , d, drop = FALSE]
  })
  
  # Call existing mat_forc for n_draws = 1, sequential
  res <- mat_forc(
    h        = h,
    n_draws  = 1L,
    n_var    = n_var,
    n_p      = n_p,
    data_    = data_,
    matrices = matrices_d,
    max_cores = 1L
  )
  
  # res$b_h: 1 x n_var x 1, res$M_h: list of length h, each n_var x n_var x 1
  b_h_vec <- drop(res$b_h[,,1])
  M_h_list <- lapply(res$M_h, function(Mj) drop(Mj[,,1]))
  
  list(b_h = b_h_vec,
       M_h = M_h_list)
}


