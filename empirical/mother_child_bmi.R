#setwd("C:/Users/alexm/NUS Dropbox/Alexander Mozdzen/MARIA-REPULSIVEWEIGHTS/2025/sip-mix/empirical")

# Load packages
library(rdist)
library(brms)
library(salso)


# Create grid
all_gammas <- c(0, 0.1, 0.5, 1, 2, 5)
all_zetas <- c(0.001, 0.1, 0.5, 1, 2, 3 )

all_gammas <- c(0, 0.5, 1, 2)
all_zetas <- c(0.001, 0.1, 0.5, 1, 2)

#all_zetas <- c(0.001, 0.1, 0.5, 1, 2, 3, 5)
all_datasets <- c("m_bmi_ch_bmi_w07", "m_bmi_ch_bmi_w10", "m_bmi_ch_bmi_w12", "m_bmi_ch_bmi_w14")

all_updates <- c("fixed")
all_alphas <- c( 1, 5) # 0. already run
repulsive_grid <- expand.grid(all_gammas, all_zetas, all_datasets, all_updates, all_alphas)
colnames(repulsive_grid) <- c("gamma", "zeta", "dataset", "update", "alpha")

dir.create(paste0("grid/"), showWarnings = FALSE, recursive = TRUE)
saveRDS(repulsive_grid, paste0("grid/","mother_child_bmi_grid.rds"))

n_save = 100;n_burn = 100;n_thin = 1


grid_eval <- function(run, repulsive_grid, n_save, n_burn, n_thin){
  
  # tryCatch({
  #setwd("C:/Users/alexander/Dropbox/MARIA-REPULSIVEWEIGHTS/Code/mixtures/fully_repulsive/repulsion_prior/SBe_prior/cluster/parallel_tryouts")
  library(rdist)
  library(brms)
  library(MASS)
  library(mclust)
  library(combinat)
  library(plyr)
  library(reshape2)
  #library(beepr)
  library(AntMAN)
  library(salso)
  library(mvtnorm)
  library(MCMCpack)
  # Source functions:
  source("src/rep_mix_main_new.R")
  source("src/log_sdir.R")
  Rcpp::sourceCpp("src/rcpp_functions.cpp")
  
  # # Assign variable based on run-id
  zetas_samp <- rep(repulsive_grid$zeta[run], 2)
  gamma_samp <- repulsive_grid$gamma[run]
  dataset <- repulsive_grid$dataset[run]
  update <- repulsive_grid$update[run]
  alpha_prior <- repulsive_grid$alpha[run]
  
  filename <- paste0( paste0(colnames(repulsive_grid[run,-3]), "_"), paste0(repulsive_grid[run,-3]), collapse = "_" )
  
  y <- as.matrix(readRDS(paste0("empirical/data/", dataset, ".rds")))
  # Set C for C for now
  C = 3
  D = 2
  N <- nrow(y)
  # Define hyperparameters and starting values:
  
  rep_prior_shape <- 0.1
  rep_prior_rate <- 1
  M_prior <- C - 1
  p_b <- 0.5
  alpha_prior <- alpha_prior
  beta_prior = alpha_prior
  sigma2_prior_shape <- 3
  sigma2_prior_rate <- 2

  # hyperparameters <- list(mu_prop_means =mu_prop_means,
  #                         mu_prop_Sigma = mu_prop_Sigma,
  #                         Sigma_prior_mat = diag(D),
  #                         Sigma_prior_df = 2,
  #                         alpha_prior = alpha_prior, beta_prior = beta_prior,
  #                         rep_prior_shape = rep_prior_shape, rep_prior_rate = rep_prior_rate,
  #                         M_prior = M_prior,
  #                         p_b = p_b,
  #                         g_to_z_ratio = 1,
  #                         g_time_z_ratio = 1)
  
  hyperparameters <- list(
    mu_prior_mean = 0, #mean(y),
    mu_prior_var_scale = 1e-05,
    #mu_prior_var_scale = 1,
    mu_prior_var = 100000,
    mu_prop_means = if(D == 1){ mean(y)} else{colMeans(y)},
    mu_prop_Sigma = if(D == 1) 1 else diag(rep(1, D)),
    Sigma_prior_mat = diag(D),
    Sigma_prior_df = 2,
    sigma2_prior_shape = sigma2_prior_shape,
    sigma2_prior_rate = sigma2_prior_rate,
    alpha_prior = alpha_prior,
    beta_prior = alpha_prior,
    M_prior = M_prior,
    p_b = 0.5,
    rep_prior_shape = rep_prior_shape,
    rep_prior_rate = rep_prior_rate,
    g_to_z_ratio = 1,
    g_time_z_ratio = 1,
    repulsive_locations = TRUE,
    mu_independent = FALSE,
    weights_prior = "selberg"
  )

  list2env(hyperparameters, envir = .GlobalEnv)

  # Define updates
  updates <- list(
    update_allocs = TRUE, update_gamma = FALSE, update_zetas = FALSE, g_to_z = FALSE, g_time_z = FALSE,
    update_var = TRUE, update_mu = TRUE, update_weight = TRUE,
    update_M = TRUE, update_post_dens = TRUE
  )
  
  list2env(updates, envir = .GlobalEnv)
  
  # Predefine M before
  if(update_M){
    M_a_samp = max(2, rpois(1, M_prior))
    M_na_samp = max(2, rpois(1, M_prior))
  }else{
    M_a_samp = C
    M_na_samp = 0
  }
  
  M_samp = M_a_samp + M_na_samp
  allocs_samp = rep(c(1:M_a_samp), length.out = N)
  mus_a_samp = t(rmvnorm(M_a_samp, mu_prop_means, mu_prop_Sigma ))
  mus_na_samp = t(rmvnorm(M_na_samp, mu_prop_means, mu_prop_Sigma ))
  mus_samp <- cbind(mus_a_samp, mus_na_samp)
  Sigmas_a_samp = array(diag(D), dim = c(D, D, M_a_samp))
  Sigmas_na_samp = array(diag(D), dim = c(D, D, M_na_samp))
  weights_samp = c(brms::rdirichlet(1, rep(alpha_prior, M_samp)))
  tuning_sd_mu = 1
  tuning_sd_gamma = 0.1
  tuning_sd_zeta = 0.1
  
  starting_values <- list(
    M_a_samp = M_a_samp,
    M_na_samp = M_na_samp,
    M_samp = M_samp,
    allocs_samp = allocs_samp,
    mus_a_samp = mus_a_samp,
    mus_na_samp = mus_na_samp,
    mus_samp = mus_samp,
    Sigmas_a_samp = Sigmas_a_samp,
    Sigmas_na_samp = Sigmas_na_samp,
    weights_samp = weights_samp,
    gamma_samp = gamma_samp,
    zetas_samp = zetas_samp,
    tuning_sd_mu = tuning_sd_mu,
    tuning_sd_gamma = tuning_sd_gamma,
    tuning_sd_zeta = tuning_sd_zeta
  )
  
  
  
  #########################################################################
  set.seed(1)
  
  # startt <- Sys.time()
  
  # results <- rep_mix_main(y, updates = updates, hyperparameters = hyperparameters , starting_values = starting_values,
  #                                   n_save = n_save, n_burn = n_burn, n_thin = n_thin)
  
  # endt <-  Sys.time() - startt 
  # endt
 
    startt <- Sys.time()
  
    results <- rep_mix_main(y, updates = updates, hyperparameters = hyperparameters , starting_values = starting_values,
                                    n_save = n_save, n_burn = n_burn, n_thin = n_thin)
  
    endt <-  Sys.time() - startt 
    endt
    
 

  # Clustering:
  
  dir.create(paste0("empirical/results_",dataset,"/allocs/"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(results$alloc_out, paste0("empirical/results_",dataset,"/allocs/","allocs_",filename,".rds"))
  
  s_out <- results$alloc_out
  a_cost <- 1
  
  # binder estimate using the salso package
  
  s_binder <- c(salso(results$alloc_out, loss=salso::binder(a = 1)))
  J_binder <- max(s_binder)
  nj_binder <- table(s_binder)
  
  dir.create(paste0("empirical/results_",dataset,"/s_binder/"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(s_binder, paste0("empirical/results_",dataset,"/s_binder/","s_binder_",filename,".rds"))
  
  dir.create(paste0("empirical/results_",dataset,"/post_dens/"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(colMeans(results$post_dens_out), paste0("empirical/results_",dataset,"/post_dens/","post_dens_",filename,".rds"))
  
  dir.create(paste0("empirical/results_",dataset,"/M_a/"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(results$M_a_out, paste0("empirical/results_",dataset,"/M_a/","M_a_",filename,".rds"))
  
  #  }, error = function(e) {
  #   message("Error at run = ", run, ": ", conditionMessage(e))
  #    list(result = NULL, error = conditionMessage(e), runtime = NA, run = run)
  # })

}

library(pbapply)
library(parallel)
parallel::detectCores()
cl <- parallel::makeCluster(parallel::detectCores())

op <- pboptions(type="timer")

# Quick run for testing
system.time(pblapply(1:nrow(repulsive_grid), grid_eval, repulsive_grid = repulsive_grid, n_save = 2e0,  n_burn = 10e0, n_thin = 2, cl = cl))

# safe_grid_eval_debug <- function(i, repulsive_grid, n_save, n_burn, n_thin) {
#   message("Running i = ", i)
#   tryCatch({
#     safe_grid_eval(i, repulsive_grid, n_save, n_burn, n_thin)
#   }, error = function(e) {
#     message("Error at i = ", i, ": ", conditionMessage(e))
#     return(NULL)
#   })
# }
#system.time(pblapply(1:nrow(repulsive_grid), safe_grid_eval_debug, repulsive_grid = repulsive_grid, n_save = 2e0,  n_burn = 10e0, n_thin = 2, cl = cl))


n_save = 5e3;  n_burn = 10e3; n_thin = 2
system.time(pblapply(1:nrow(repulsive_grid), grid_eval, repulsive_grid = repulsive_grid, n_save = n_save,  n_burn = n_burn, n_thin = n_thin, cl = cl))

parallel::stopCluster(cl)

#############################################################################################
