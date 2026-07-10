# detach("APRScenario", unload = TRUE, character.only = TRUE)
remotes::install_github("titoqu/APRScenario")

devtools::load_all()


# Three variables example testing running time to compute big_b and big_M (forecast matrices) compared to whole scenarios function.
# This is to assess if big_b_and_M is the bottleneck of the function or if also the next part of the scenarios code should be optimized
{
  data(NKdata)
  
  #params
  S=100 # Number of draws
  h=10 #forecast horizon
  p=10 #lags in the VAR model
  
  spec <- bsvarSIGNs::specify_bsvarSIGN$new(as.matrix(NKdata[,2:4]), p = p)
  posterior <- bsvars::estimate(spec, S = S)
  matrices<-gen_mats(posterior = posterior, specification = spec)
  
  start_time_bM <- Sys.time()
  tmp <- big_b_and_M(
    h = h,
    n_draws = dim(posterior$posterior$B)[3],
    n_var = dim(posterior$posterior$B)[1],
    n_p = (dim(posterior$posterior$A)[2] - 1) / dim(posterior$posterior$B)[1],
    matrices = matrices,
    # n_cores = 1
  )
  end_time_bM<- Sys.time()
  elapsed_time_bM<- end_time_bM- start_time_bM
  
  
  
  start_time_scen <- Sys.time()
  # and having posterior object
  scenario_result <- scenarios(h = h, 
                               path = rep(1,h), 
                               obs = 1, 
                               free_shocks = NA, 
                               posterior = posterior, 
                               matrices = matrices)
  end_time_scen <- Sys.time()
  elapsed_time_scen<- end_time_scen- start_time_scen
  
  paste0("Run time(old code) - big_b and big_M: ",sprintf("%.2f", elapsed_time_bM)," sec. - Whole scenario: ",
         sprintf("%.2f", elapsed_time_scen), " (",
         sprintf("%.2f", as.numeric(elapsed_time_scen)/as.numeric(elapsed_time_bM))," times SLOWER)"
  )
  
}
# with different set of parameters the answer is similar, the time for computing is extremely close hence it is a clear bottleneck



# Three variables example testing running time to compute b and big_M (forecast matrices) with new parallelization method. 
# This is the computational bottleneck of the scenarios() fun.)
{
  data(NKdata)
  
  #params
  S=1000 # Number of draws
  h=24 #forecast horizon
  p=10 #lags in the VAR model
  
  spec <- bsvarSIGNs::specify_bsvarSIGN$new(as.matrix(NKdata[,2:4]), p = p)
  posterior <- bsvars::estimate(spec, S = S)
  matrices<-gen_mats(posterior = posterior, specification = spec)
  
  # old function
  start_time_old <- Sys.time()
  tmp <- big_b_and_M(
    h = h,
    n_draws = dim(posterior$posterior$B)[3],
    n_var = dim(posterior$posterior$B)[1],
    n_p = (dim(posterior$posterior$A)[2] - 1) / dim(posterior$posterior$B)[1],
    matrices = matrices,
    # n_cores = 1
  )
  end_time_old <- Sys.time()
  elapsed_time_old <- end_time_old - start_time_old
  
  # new function with one core (different structure hence is should be already faster)
  start_time1 <- Sys.time()
  tmp <- big_b_and_M_pl(
    h = h,
    n_draws = dim(posterior$posterior$B)[3],
    n_var = dim(posterior$posterior$B)[1],
    n_p = (dim(posterior$posterior$A)[2] - 1) / dim(posterior$posterior$B)[1],
    matrices = matrices,
    n_cores = 1
  )
  end_time1 <- Sys.time()
  elapsed_time_1<- end_time1- start_time1
  
  # new function with max cores (8 in my machine) (should be the fastest option at least for large n_draws and lags p)
  start_time <- Sys.time()
  tmp <- big_b_and_M_pl(
    h = h,
    n_draws = dim(posterior$posterior$B)[3],
    n_var = dim(posterior$posterior$B)[1],
    n_p = (dim(posterior$posterior$A)[2] - 1) / dim(posterior$posterior$B)[1],
    matrices = matrices,
    n_cores = 8
  )
  end_time <- Sys.time()
  elapsed_time_8 <- end_time - start_time
  
  paste0("Run time - old code: ",sprintf("%.2f", elapsed_time_old)," sec. - New code with 1 core: ",
         sprintf("%.2f", elapsed_time_1)," sec. - New code with 8 cores: ",sprintf("%.2f", elapsed_time_8)," sec. (",
         sprintf("%.1f", as.numeric(elapsed_time_1)/as.numeric(elapsed_time_8))," times faster than 1 core and ",
         sprintf("%.1f", as.numeric(elapsed_time_old)/as.numeric(elapsed_time_8)), " faster than the old code.)"
         )
}


# Compare Old and New scenarios() fun. - (on Windows) the parallelized function (New) outperforms heavily the old one 
# in terms of running time (even by more than a factor of n_cores, due to new design). Especially when increasing complexity or number of draws
{
  data(NKdata)
  
  #params
  S=100 # Number of draws
  h=20 #forecast horizon
  p=10 #lags in the VAR model
  
  spec <- bsvarSIGNs::specify_bsvarSIGN$new(as.matrix(NKdata[,2:4]), p = p)
  posterior <- bsvars::estimate(spec, S = S)
  matrices<-gen_mats(posterior = posterior, specification = spec)
  
  start_time_old <- Sys.time()
  scenario_result_old <- scenarios(h = h, 
                               path = rep(1,h), 
                               obs = 1, 
                               free_shocks = NA, 
                               posterior = posterior, 
                               matrices = matrices)
  end_time_old <- Sys.time()
  elapsed_time_old<- end_time_old- start_time_old
  
  
  
  start_time_new <- Sys.time()
  
  scenario_result_new <- scenarios_pl(h = h, 
                               path = rep(1,h), 
                               obs = 1, 
                               free_shocks = NA, 
                               posterior = posterior, 
                               matrices = matrices,
                               n_cores = 8)
  end_time_new <- Sys.time()
  elapsed_time_new<- end_time_new- start_time_new
  
  paste0("Run time (scenarios) - Old code: ",sprintf("%.2f", elapsed_time_old)," sec. - New code: ",
         sprintf("%.2f", elapsed_time_new), " (",
         sprintf("%.1f", as.numeric(elapsed_time_old)/as.numeric(elapsed_time_new))," times FASTER)"
  )
  
}


#Check if there are differences in the two scenarios() results (no differences found at 12 decimal points precision)
{
  old <- unlist(scenario_result_old, recursive = TRUE, use.names = FALSE)
  new <- unlist(scenario_result_new, recursive = TRUE, use.names = FALSE)
  
  different <- round(old, 12) != round(new, 12)
  
  paste0(
    "Different at 12 decimals: ",
    sum(different, na.rm = TRUE), " / ", length(different),
    " (", sprintf("%.2f", 100 * mean(different, na.rm = TRUE)), "%)",
    "\nMaximum absolute difference: ",
    sprintf("%.12f", max(abs(old - new), na.rm = TRUE))
  )
}

