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
  
  K_0 <- diag(n_var)
  K_h <- vector("list", h)
  K_h[[1]] <- K_0
  
  if (h > 1) {
    for (i in 2:h) {
      tmp2 <- matrix(0, n_var, n_var)
      
      for (j in 1:(i - 1)) {
        if (j <= n_p) {
          tmp1 <- matrices$B_list[[j]][, , d]
        } else {
          tmp1 <- matrix(0, n_var, n_var)
        }
        
        tmp2 <- tmp2 + K_h[[i - j]] %*% tmp1
      }
      
      K_h[[i]] <- tmp2 + K_0
    }
  }
  
  M_h <- vector("list", h)
  M_h[[1]] <- matrices$M[, , d]
  
  if (h > 1) {
    for (i in 2:h) {
      tmp2 <- matrix(0, n_var, n_var)
      
      for (j in 1:min(i - 1, n_p)) {
        tmp2 <- tmp2 + M_h[[j]] %*% matrices$B_list[[j]][, , d]
      }
      
      M_h[[i]] <- tmp2
    }
  }
  
  N_p_list <- vector("list", n_p)
  
  for (l in seq_len(n_p)) {
    tmp00 <- vector("list", h)
    tmp00[[1]] <- matrices$B_list[[l]][, , d]
    
    if (h > 1) {
      for (i in 2:h) {
        tmp2 <- matrix(0, n_var, n_var)
        
        for (j in 1:min(i - 1, n_p)) {
          tmp2 <- tmp2 +
            tmp00[[i - j]] %*%
            matrices$B_list[[j]][, , d]
        }
        
        if ((l + i - 1) <= n_p) {
          tmp <- matrices$B_list[[l + i - 1]][, , d]
        } else {
          tmp <- matrix(0, n_var, n_var)
        }
        
        tmp00[[i]] <- tmp2 + tmp
      }
    }
    
    N_p_list[[l]] <- tmp00
  }
  
  b_all <- matrix(0, n_var, h)
  
  for (hh in seq_len(h)) {
    b_hh <- as.numeric(matrices$intercept[, d] %*% K_h[[hh]])
    
    for (cnt in seq_len(n_p)) {
      y_lag <- data_[
        (1 + n_var * (cnt - 1)):(n_var * cnt),
        ncol(data_)
      ]
      
      b_hh <- b_hh +
        as.numeric(t(y_lag) %*% N_p_list[[cnt]][[hh]])
    }
    
    b_all[, hh] <- b_hh
  }
  
  list(
    b_h = as.vector(b_all),
    M_h = M_h
  )
}