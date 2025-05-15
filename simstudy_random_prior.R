#setwd("~/paper/fullrep/M_random")
#setwd("~/paper/fullrep/M_random")

#rm(list = ls())

#Create grid:
#all_alphas <- c(0.05, 0.1, 1, 10)
#all_alphas <- c(0.1, 1)
all_alphas <- c( 1)
#all_gammas <- c(0, 0.25, 0.5,  1, 98, 100, 101)
#all_zetas <- c(0.001,  1,  3, 5, 98)

# pre-2025:
all_gammas <- c(100)
all_zetas <- c(0.001,  1, 5, 100)

# 2025:
all_gammas <- c(0, 0.25, 1, 98, 100)
all_zetas <- c(0.1, 99)

repulsive_grid <- expand.grid(all_alphas, all_gammas, all_zetas)
colnames(repulsive_grid) <- c("alpha", "gamma", "zeta")

dir.create(paste0("simstudy/grid/"), showWarnings = FALSE, recursive = TRUE)
saveRDS(repulsive_grid, paste0("simstudy/grid/","repulsive_grid.rds"))

if(FALSE){
  run = 1
  dataset = "simdata_1"
  n_save = 5e2;  n_burn = 10e2; n_thin = 1
}

grid_eval <- function(run, repulsive_grid, n_save, n_burn, n_thin, dataset){
  
  # Load libraries
  library(rdist)
  library(brms)
  library(MASS)
  library(mclust)
  library(combinat)
  library(plyr)
  library(reshape2)
  library(salso)
  library(mvtnorm)
  library(MCMCpack)
  # Source functions
  source("src/log_sdir.R")
  source("src/rep_mix_main_revised.R")
  # 2 has the option of an independent prior for mu
  Rcpp::sourceCpp("src/rcpp_functions.cpp")
  
  filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
  
  if(dataset == "simdata_1"){
    
    set.seed(3)
    N <- 300
    D = 2
    true_weights <- c(0.2, 0.2, 0.2, 0.3, 0.1)
    true_mus_1 <- c(-3, -3, 3, 3, -1)
    true_mus_2 <- c(-2.5, 3, -3, 3, 0)
    true_mus <- rbind(true_mus_1, true_mus_2)
    C <- length(true_mus_1)
    
    true_Sigmas <- array(
      matrix(c(3, 1, 1, 3), 2, 2),
      dim = c(2, 2, C)
    )
    # elongate 4th cluster to overlap with the 3rd:
    #true_Sigmas[,,4] <- matrix(c(1, 0, 0, 10), 2, 2)
    # make redundant cluster smaller:
    true_Sigmas[,,5] <- matrix(c(0.25, 0, 0, 0.25), 2, 2)
    
    true_allocs <- sample(C, N, prob = true_weights, replace = TRUE)
    
    y <- matrix(NA, N, D)
    for(i in 1:N){
      
      y[i,] <- rmvnorm(1, mean = true_mus[,true_allocs[i]], sigma = true_Sigmas[,,true_allocs[i]])
      
    }
    
    plot(y)
    #points(y, rep(0, N))
    
    dir.create(paste0("simstudy/datasets/",dataset), showWarnings = FALSE, recursive = TRUE)
    saveRDS(y, paste0("simstudy/datasets/",dataset,".rds"))
    
  }else  if(dataset == "simdata_2"){
    
    
    
  }
  
  # # Assign variable based on run-id
  zetas_samp <- rep(repulsive_grid$zeta[run], D)
  gamma_samp <- repulsive_grid$gamma[run]
  alpha_prior <- repulsive_grid$alpha[run]
  
  sigma2_prior_shape <- 3
  sigma2_prior_rate <- 2
  #sigma2_prior_shape <- 2.11
  #sigma2_prior_rate <- 0.11
  rep_prior_shape <- 3
  rep_prior_rate <- 3
  M_prior <- C - 1
  p_b <- 0.5
  
  # We can define this globally since sometimes its just turned off through g_to_z
  g_to_z_ratio = zetas_samp[1]
  
  if(gamma_samp == 98){
    gamma_samp = 1
    update_gamma = TRUE
    rep_prior_shape = 3
    rep_prior_rate = 2
    g_to_z = FALSE
  }else if(gamma_samp == 99){
    gamma_samp = 1
    update_gamma = TRUE
    rep_prior_shape = 12
    rep_prior_rate = 8
    g_to_z = FALSE
  }else if(gamma_samp == 100){
    gamma_samp = 1
    update_gamma = TRUE
    rep_prior_shape = 3
    rep_prior_rate = 2
    g_to_z = TRUE
  }else if(gamma_samp == 101){
    gamma_samp = 1
    update_gamma = FALSE
    rep_prior_shape = 12
    rep_prior_rate = 8
    g_to_z = TRUE
  }else{
    update_gamma = FALSE
    g_to_z = FALSE
  }  
  
  if(zetas_samp[1] == 98){
    zeta_samp = rep(1, D)
    update_zetas = TRUE
    rep_prior_shape = 3
    rep_prior_rate = 2
  }else if(zetas_samp[1] == 99){
    zeta_samp = rep(1, D)
    update_zetas = TRUE
    rep_prior_shape = 0.1
    rep_prior_rate = 1
  }else{
    update_zetas = FALSE
  }
  
  hyperparameters <- list(
    mu_prior_mean = 0, #mean(y),
    mu_prior_var_scale = 1e-05,
    #mu_prior_var_scale = 1,
    mu_prior_var = 100000,
    mu_prop_means = if(D == 1){ mean(y)} else{colMeans(y)},
    mu_prop_Sigma = if(D == 1) 1 else diag(rep(10, D)),
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
    g_to_z_ratio = g_to_z_ratio,
    g_time_z_ratio = 10,
    repulsive_locations = TRUE,
    mu_independent = FALSE,
    weights_prior = "selberg"
  )
  
  list2env(hyperparameters, envir = .GlobalEnv)
  
  updates <- list(
    update_allocs = TRUE, update_gamma = update_gamma, update_zetas = update_zetas, 
    g_to_z = g_to_z, g_time_z = FALSE,
    update_var = TRUE, update_mu = TRUE, update_weight = TRUE,
    update_M = TRUE, update_post_dens = TRUE
  )
  
  list2env(updates, envir = .GlobalEnv)
  
  ##############################################################################
  
  if(update_allocs & update_M){
    M_a_samp = C
    M_na_samp = max(2, rpois(1, M_prior))
    M_samp = M_a_samp + M_na_samp
    allocs_samp = rep(c(1:M_a_samp), length.out = N)
  }else if(!update_allocs){
    allocs_samp = true_allocs
    M_a_samp = C
    M_na_samp = 0
    M_samp = M_a_samp 
  }else if(update_allocs & !update_M){
    M_a_samp = C
    allocs_samp = rep(c(1:M_a_samp), length.out = N)
    M_na_samp = 0
    M_samp = M_a_samp
  }
  
  if(update_mu){
    if(D == 1){
      mus_a_samp = rnorm(M_a_samp, mu_prop_means, mu_prop_Sigma)
      mus_na_samp = rnorm(M_na_samp, mu_prop_means, mu_prop_Sigma)
      mus_samp <-  c(mus_a_samp, mus_na_samp)
    }else{
      mus_na_samp = t(rmvnorm(M_na_samp, mu_prop_means, mu_prop_Sigma))
      mus_a_samp = t(rmvnorm(M_a_samp, mu_prop_means, mu_prop_Sigma))
      mus_samp <- cbind(mus_a_samp, mus_na_samp)
    }
  }else{
    mus_a_samp = true_mus
    mus_na_samp = rep(0, M_na_samp)
    mus_samp <- c(mus_a_samp, mus_na_samp)
  }
  
  if(update_weight){
    weights_samp = c(brms::rdirichlet(1, rep(alpha_prior, M_samp)))
  }else{
    weights_samp = c(true_weights, rep(0, M_na_samp))
  }
  
  if(update_var){
    if(D == 1){
      Sigmas_a_samp = rep(1, M_a_samp)
      Sigmas_na_samp = rep(1, M_na_samp)
      Sigmas_samp = c(Sigmas_a_samp, Sigmas_na_samp)
    }else{
      Sigmas_a_samp = array(diag(D), dim = c(D, D, M_a_samp))
      Sigmas_na_samp = array(diag(D), dim = c(D, D, M_na_samp))
      Sigmas_samp <- array(c(Sigmas_a_samp, Sigmas_na_samp), dim = c(D, D, M_samp))
    }
  }else{
    Sigmas_a_samp = true_sds^2
    Sigmas_na_samp = numeric(0)
    Sigmas_samp = Sigmas_a_samp
  }
  
  
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
  
  true_weights <- c(0.2, 0.2, 0.2, 0.3, 0.1)
  
  #########################################################################
  set.seed(1)
  
  startt <- Sys.time()
  
  results <- rep_mix_main(y, updates = updates, hyperparameters , starting_values = starting_values,
                          n_save = n_save, n_burn = n_burn, n_thin = n_thin, true_weights)
  
  endt <-  Sys.time() - startt 
  endt
  
  # Quick Check:
  # n_save = 5e3;  n_burn = 5e3; n_thin = 1
  # (results$mu_out %>% 
  #   data.frame() %>% 
  #   dplyr::mutate(iter = 1:n_save) %>% 
  #   pivot_longer(cols = -iter, names_to = "component", values_to = "value") %>% 
  #   ggplot() + 
  #   geom_density(alpha = 0.4, aes(x = value, col = component, fill = component),stat = "density") +
  #   geom_vline(xintercept = true_mus)) / (results$mu_out %>% 
  #   data.frame() %>% 
  #   dplyr::mutate(iter = 1:n_save) %>% 
  #   pivot_longer(cols = -iter, names_to = "component", values_to = "value") %>% 
  #   ggplot(aes(x = iter, y = value,  col = component, fill = component)) + 
  #   geom_line() +
  #   geom_hline(yintercept = true_mus))
  # 
  # results$weights_out %>% 
  #   data.frame() %>% 
  #   dplyr::mutate(iter = 1:n_save) %>% 
  #   pivot_longer(cols = -iter, names_to = "component", values_to = "value") %>% 
  #   ggplot() + 
  #   geom_density(alpha = 0.4, aes(x = value, col = component, fill = component),stat = "density") 
  # 
  # results$weights_out %>%
  #   data.frame() %>%
  #   dplyr::mutate(iter = 1:n_save) %>%
  #   pivot_longer(cols = -iter, names_to = "component", values_to = "value") %>%
  #   ggplot(aes(x = iter, y = value,  col = component, fill = component)) +
  #   geom_line()
  # 
  # results$Sigma_out %>% 
  #   data.frame() %>% 
  #   dplyr::mutate(iter = 1:n_save) %>% 
  #   pivot_longer(cols = -iter, names_to = "component", values_to = "value") %>% 
  #   ggplot() + 
  #   geom_density(alpha = 0.4, aes(x = value, col = component, fill = component),stat = "density") +
  #   geom_vline(xintercept = true_sds^2)
  # 
  # results$Sigma_out %>% 
  #   data.frame() %>% 
  #   dplyr::mutate(iter = 1:n_save) %>% 
  #   pivot_longer(cols = -iter, names_to = "component", values_to = "value") %>% 
  #   ggplot(aes(x = iter, y = value,  col = component, fill = component)) + 
  #   geom_line() +
  #   geom_hline(yintercept = true_sds^2)
  
  # post dens
  
  dir.create(paste0("simstudy/results/",dataset,"/post_dens/"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(colMeans(results$post_dens_out), paste0("simstudy/results/",dataset,"/post_dens/","post_dens_",filename,".rds"))
  
  # M_a
  
  dir.create(paste0("simstudy/results/",dataset,"/M_a/"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(results$M_a_out, paste0("simstudy/results/",dataset,"/M_a/","M_a_",filename,".rds"))
  
  # allocations:
  
  dir.create(paste0("simstudy/results/",dataset,"/allocs/"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(results$alloc_out, paste0("simstudy/results/",dataset,"/allocs/","allocs_",filename,".rds"))
  
  # Mus
  
  dir.create(paste0("simstudy/results/",dataset,"/mus/"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(results$mu_out, paste0("simstudy/results/",dataset,"/mus/","mus_",filename,".rds"))
  
  # weights
  
  dir.create(paste0("simstudy/results/",dataset,"/weights/"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(results$weights_out, paste0("simstudy/results/",dataset,"/weights/","weights_",filename,".rds"))
  
  # sigmas
  
  dir.create(paste0("simstudy/results/",dataset,"/sigmas/"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(results$Sigma_out, paste0("simstudy/results/",dataset,"/sigmas/","sigmas_",filename,".rds"))
  
  # gammas
  
  dir.create(paste0("simstudy/results/",dataset,"/gammas/"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(results$gamma_out, paste0("simstudy/results/",dataset,"/gammas/","gammas_",filename,".rds"))
  
  # zetas
  
  dir.create(paste0("simstudy/results/",dataset,"/zetas/"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(results$zetas_out, paste0("simstudy/results/",dataset,"/zetas/","zetas_",filename,".rds"))
  
  # s_out <- results$alloc_out
  # a_cost <- 1
  # 
  # # binder estimate using the salso package
  # 
  # s_binder <- c(salso(results$alloc_out, loss=salso::binder(a = 1)))
  # J_binder <- max(s_binder)
  # nj_binder <- table(s_binder)
  # 
  # dir.create(paste0("simstudy/results_both_",dataset,"/s_binder/"), showWarnings = FALSE, recursive = TRUE)
  # saveRDS(s_binder, paste0("simstudy/results_both_",dataset,"/s_binder/","s_binder_",filename,".rds"))
  # 
  # dir.create(paste0("simstudy/results_both_",dataset,"/post_dens/"), showWarnings = FALSE, recursive = TRUE)
  # saveRDS(colMeans(results$post_dens_out), paste0("simstudy/results_both_",dataset,"/post_dens/","post_dens_",filename,".rds"))
  # 
  # dir.create(paste0("simstudy/results_both_",dataset,"/M_a/"), showWarnings = FALSE, recursive = TRUE)
  # saveRDS(results$M_a_out, paste0("simstudy/results_both_",dataset,"/M_a/","M_a_",filename,".rds"))
  
}

library(pbapply)
library(parallel)
parallel::detectCores()
cl <- parallel::makeCluster(parallel::detectCores())

op <- pboptions(type="timer")

# Quick run for testing
#system.time(pblapply(1:nrow(repulsive_grid), grid_eval, repulsive_grid = repulsive_grid, n_save = 2e1,  n_burn = 10e1, n_thin = 2, cl = cl))

n_save = 5e3;  n_burn = 5e3; n_thin = 1

system.time(pblapply(1:nrow(repulsive_grid), grid_eval, repulsive_grid = repulsive_grid, n_save = n_save,  n_burn = n_burn, n_thin = n_thin, cl = cl, dataset = "simdata_1"))
#system.time(pblapply(1, grid_eval, repulsive_grid = repulsive_grid, n_save = n_save,  n_burn = n_burn, n_thin = n_thin, cl = cl, dataset = "simdata_1"))
#system.time(pblapply(which(repulsive_grid$gamma == 101), grid_eval, repulsive_grid = repulsive_grid, n_save = n_save,  n_burn = n_burn, n_thin = n_thin, cl = cl, dataset = "simdata_1"))

parallel::stopCluster(cl)


























plots = FALSE

# variance_cluster_sizes <- function(cluster_matrix) {
#   apply(cluster_matrix, 1, function(clusters) {
#     sizes <- as.numeric(table(clusters))
#     var(sizes)
#   })
# }
# 
# gini_coefficient <- function(cluster_matrix) {
#   apply(cluster_matrix, 1, function(clusters) {
#     sizes <- table(clusters)
#     ineq(sizes, type = "Gini")
#   })
# }
# 
# entropy <- function(cluster_matrix) {
#   apply(cluster_matrix, 1, function(clusters) {
#     sizes <- table(clusters)
#     p <- sizes / sum(sizes)
#     -sum(p * log(p), na.rm = TRUE)
#   })
# }
# 
# 
#setwd("~/paper/fullrep/M_random")
dataset <- "simdata_1"
y <- readRDS(paste0("simstudy/datasets/",dataset,".rds"))
repulsive_grid <- readRDS(paste0("simstudy/grid/repulsive_grid.rds"))

###############################################################################
# Plots
###############################################################################

from_to <- apply(y, 2, range)
y_line1 <- seq(from_to[1,1], from_to[2,1], length.out = 2e1)
y_line2 <- seq(from_to[1,2], from_to[2,2], length.out = 2e1)
y_grid <- as.matrix(expand.grid(y_line1, y_line2))

#################################################################################
# For convenience
#################################################################################
dir.create(paste0("simstudy/results/",dataset,"/plots/"), showWarnings = FALSE, recursive = TRUE)
#library(label.switching)
library(mcclust.ext)
library(ggalt)
all_mus.w <- data.frame()
all_weights.w <- data.frame()
all_sigmas.w <- data.frame()
all_M_a.w <- data.frame()
all_post_dens.w <- data.frame()
all_psms <- list()
binder.ggs <- list()
vi.ggs <- list()

all_ws.l_list <- list()


#detach("package:plyr", unload=TRUE)

for(run in 1:nrow(repulsive_grid)){
#for(run in which(repulsive_grid$gamma == 101)){
  print(run)
  filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )

  w3.w <- matrix(NA, n_save, 3)
  for(i in 1:n_save){
    w3.w[i,] <- readRDS(paste0("simstudy/results/",dataset,"/weights/","weights_",filename,".rds"))[[i]][1:3]
  }

  # weights:
  all_ws.l_list[[run]] <-
    cbind(repulsive_grid[run,], w3.w) %>% mutate(iter = 1:n_save) %>%
    pivot_longer(cols = -c(colnames(repulsive_grid),"iter"), names_to = "component", values_to = "value")

  # M_a

  M_a.w <- cbind(repulsive_grid[run,],
                 M_a = readRDS(paste0("simstudy/results/",dataset,"/M_a/","M_a_",filename,".rds")))%>%
                    mutate(iter = 1:n_save)

  all_M_a.w <- rbind(all_M_a.w, M_a.w)

  zetas_samp <- rep(repulsive_grid$zeta[run], 2)
  gamma_samp <- repulsive_grid$gamma[run]
  alpha_prior <- repulsive_grid$alpha[run]

  allocs <- readRDS(paste0("simstudy/results/",dataset,"/allocs/","allocs_",filename,".rds"))
  s_binder <- c(salso(allocs, loss=salso::binder(a = 1)))

  df <- data.frame(x = y[,1], y = y[,2], Cluster = factor(s_binder))

  binder.ggs[[run]] <-
    ggplot() +
    geom_point(data = df, aes(x = x, y = y, color = Cluster),  size = 2) +
    geom_encircle(data = df, aes(x = x, y = y, group = Cluster, fill = Cluster),
                  size = 0.5, expand = 0, alpha = 0.15, s_shape = 0.5,  spread=0.025) +
    labs(title = paste0("Binder \n zeta = ", zetas_samp[1], " and gamma = ", gamma_samp), x = "BMI", y = "CEBQ Person Parameter", legend = "test") +
    theme_minimal()

  s_vi <- c(salso(allocs, loss=salso::VI(a = 1)))

  df <- data.frame(x = y[,1], y = y[,2], Cluster = factor(s_vi))

  vi.ggs[[run]] <-
    ggplot() +
    geom_point(data = df, aes(x = x, y = y, color = Cluster),  size = 2) +
    geom_encircle(data = df, aes(x = x, y = y, group = Cluster, fill = Cluster),
                  size = 0.5, expand = 0, alpha = 0.15, s_shape=0.5,  spread=0.025) +
    labs(title = paste0("vi \n zeta = ", zetas_samp[1], " and gamma = ", gamma_samp), x = "BMI", y = "CEBQ Person Parameter", legend = "test") +
    theme_minimal()


}

all_ws.l <- do.call(rbind, all_ws.l_list)
# Double checking which runs had good mixing:
ess_threshold <- 50
#
all_ws_filtered <- all_ws.l %>%
  na.omit() %>%
  group_by(alpha, gamma, zeta, component) %>%
  summarise(ess = effectiveSize(value), .groups = 'drop')
#
# constellations_majority_ess <- all_ws_filtered %>%
#   group_by(alpha, gamma,zeta) %>%
#   summarise(num_low_ess_components = sum(ess < ess_threshold),  .groups = 'drop') %>%
#   # filter(num_low_ess_components < ceiling(M_samp / 2))
#   filter(num_low_ess_components < 2)
#
# ess_ws.l <- all_ws.l %>%
#   inner_join(constellations_majority_ess, by = c("alpha", "gamma", "zeta"))
#
# nrow(ess_ws.l) / nrow(all_ws.l)
#
# ess_ws.l %>%
#   filter(alpha ==1) %>%
#   ggplot(aes(x = iter, y = value, col = component))+
#   geom_line() +
#   facet_wrap(zeta ~  gamma, ncol = length(all_gammas))

# all_ws.l%>%
#   filter(alpha ==1) %>%
#   ggplot(aes(x = iter, y = value, col = component))+
#   geom_line() +
#   facet_wrap(zeta ~  gamma, ncol = length(all_gammas))

# ---> all look decent!


# Scatterplot of the data with true clusters

# set.seed(3)
# N <- 300
# D = 2
# true_weights <- c(0.2, 0.2, 0.2, 0.3, 0.1)
# true_mus_1 <- c(-3, -3, 3, 3, -1)
# true_mus_2 <- c(-2.5, 3, -3, 3, 0)
# true_mus <- rbind(true_mus_1, true_mus_2)
# C <- length(true_mus_1)
#
# true_Sigmas <- array(
#   matrix(c(3, 1, 1, 3), 2, 2),
#   dim = c(2, 2, C)
# )
# # elongate 4th cluster to overlap with the 3rd:
# #true_Sigmas[,,4] <- matrix(c(1, 0, 0, 10), 2, 2)
# # make redundant cluster smaller:
# true_Sigmas[,,5] <- matrix(c(0.25, 0, 0, 0.25), 2, 2)
set.seed(3)
true_allocs <- sample(C, N, prob = true_weights, replace = TRUE)
y2 <- data.frame(y, cluster = true_allocs )

# Scatter plot with clusters encircled


#pdf(paste0("simstudy/results/",dataset,"/plots/data_cluster.pdf"), width = 10, height = 10)

ggplot(y2, aes(x = X1, y = X2, color = factor(cluster))) +
  geom_point(size = 3) +
  geom_encircle(aes(group =factor( cluster),  fill = factor(cluster)),
                size = 0.5, expand = 0, alpha = 0.15, s_shape = 0.5,  spread=0.025)+
  labs( fill = "Cluster",
        color = "Cluster",
        x = "",
        y = "") +
  theme_minimal() +
  scale_fill_manual(values = c(pal_npg()(10)[-2], pal_npg()(10)[2])) +
  scale_color_manual(values = c(pal_npg()(10)[-2], pal_npg()(10)[2]))

#dev.off()

# "PSM" of true clustering:
library(ggplot2)
library(reshape2)

# Step 1: Create the true similarity matrix
N <- length(true_allocs)
true_similarity <- matrix(0, nrow = N, ncol = N)

for (i in 1:N) {
  for (j in 1:N) {
    true_similarity[i, j] <- ifelse(true_allocs[i] == true_allocs[j], 1, 0)
  }
}

# Step 2: Convert to data frame
true_similarity_df <- melt(true_similarity)
colnames(true_similarity_df) <- c("Observation1", "Observation2", "Similarity")

# Perform hierarchical clustering on the similarity matrix
hc <- hclust(as.dist(1 - true_similarity))
reordered_indices <- hc$order
# true_similarity_df$SameCluster <- ifelse(true_allocs[true_similarity_df$Observation1] ==
#                                            true_allocs[true_similarity_df$Observation2],
#                                          true_allocs[true_similarity_df$Observation1], NA)


true_similarity_df$Observation1 <- factor(true_similarity_df$Observation1, levels = reordered_indices)
true_similarity_df$Observation2 <- factor(true_similarity_df$Observation2, levels = reordered_indices)

head(true_similarity_df)

# Calculate cluster centers
cluster_centers <- data.frame(Observation = 1:length(true_allocs), Cluster = true_allocs)
cluster_centers <- aggregate(Observation ~ Cluster, data = cluster_centers, FUN = median)

cluster_centers$Observation[1] <- 10
cluster_centers$Observation[2] <- 50
cluster_centers$Observation[3] <- 110
cluster_centers$Observation[4] <- 190
cluster_centers$Observation[5] <- 270

cluster_centers$Cluster <- c(5,1,3,4,2)

#true_similarity_df$Cluster <- as.factor(true_allocs[true_similarity_df$Observation1])

pdf(paste0("simstudy/results/",dataset,"/plots/true_similarity_matrix.pdf"), width = 10, height = 10)

ggplot(true_similarity_df, aes(x = Observation1, y = Observation2, fill = Similarity)) +
  geom_tile() +
  scale_fill_manual(values = pal_npg()(length(unique(true_allocs)))) +
  #scale_fill_manual(values = c(pal_npg()(10)[-2], pal_npg()(10)[2])) +
  #scale_color_manual(values = c(pal_npg()(10)[-2], pal_npg()(10)[2]))
  scale_fill_gradient(low = "#f4f3ff", high = pal_npg()(10)[7]) +
  labs(title = "", x = "", y = "", fill = "") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  geom_text(data = cluster_centers, aes(x = Observation, y = Observation, label = Cluster), inherit.aes = FALSE, color = "black", size = 10)

dev.off()



######################################################################################
# PSM ggplot
# 4x4 psms done with ggplot:
######################################################################################

# Select those psms that are needed manually. Currently:
# zeta: 0.001, 0.1, 5, G(0.1,1) = 99
# gamma: 0, 0.25, 1, G(3, 2) = 98

repulsive_grid

repulsive_grid_fixed <- data.frame(repulsive_grid) %>% filter(!(gamma %in% c(0.5,99:101)) & !(zeta %in% c(98, 1, 3))) %>% arrange(gamma, zeta)

gg_psms <- list()
for(run in 1:nrow(repulsive_grid_fixed)){

  #filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
  filename <- paste0( paste0(colnames(repulsive_grid_fixed[run,]), "_"), paste0(repulsive_grid_fixed[run,]), collapse = "_" )
  allocs <- readRDS(paste0("simstudy/results/",dataset,"/allocs/allocs_",filename,".rds"))

  #ind <- which(run == psm_fixed_plot_order)
  ind <- run

  new_color <- rep(1:4, each = 4)

  # Compute the Posterior Similarity Matrix from the allocation matrix
  posterior_similarity <- matrix(0, nrow = ncol(allocs), ncol = ncol(allocs))

  for (i in 1:nrow(allocs)) {
    posterior_similarity <- posterior_similarity + (outer(allocs[i, ], allocs[i, ], "=="))
  }

  # Normalize the similarity matrix by the number of iterations
  posterior_similarity <- posterior_similarity / nrow(allocs)

  # Convert the matrix to a format suitable for ggplot2
  library(reshape2)
  similarity_df <- melt(posterior_similarity)
  colnames(similarity_df) <- c("Observation1", "Observation2", "Similarity")

  # Perform hierarchical clustering on the similarity matrix
  hc <- hclust(as.dist(1 - posterior_similarity))
  reordered_indices <- hc$order
  similarity_df$Observation1 <- factor(similarity_df$Observation1, levels = reordered_indices)
  similarity_df$Observation2 <- factor(similarity_df$Observation2, levels = reordered_indices)

  # Very Light Gray: #e0e0e0
  # Ivory: #fffff0
  # Lavender: #f4f3ff
  # Beige: #f5f5dc
  # Mint Cream: #f5fff5


  # good: "#f0f0f0"

  # Plot the clustered heatmap
  gg_psms[[ind]] <- ggplot(similarity_df, aes(x = Observation1, y = Observation2, fill = Similarity)) +
    geom_tile() +
    scale_fill_gradient(low = "#f4f3ff", high = pal_npg()(4)[new_color[ind]]) +
    #scale_fill_gradient(low = "white", high = "darkred") +
    labs(title = "",
         x = "",
         y = "",
         fill = "") +
    theme_minimal() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

}


par(mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))

grid_plot <-
((gg_psms[[1]] | gg_psms[[2]] | gg_psms[[3]] | gg_psms[[4]]) &
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))/
 ( (gg_psms[[5]] | gg_psms[[6]] | gg_psms[[7]] | gg_psms[[8]]) &
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))/
  ((gg_psms[[9]] | gg_psms[[10]] | gg_psms[[11]] | gg_psms[[12]]) &
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))/
  ((gg_psms[[13]] | gg_psms[[14]] | gg_psms[[15]] | gg_psms[[16]]) &
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))) +
  plot_layout(guides = "collect") &
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

grid_plot

#pdf(paste0("simstudy/results/",dataset,"/plots/true_similarity_matrix.pdf"), width = 10, height = 10)
png( paste0("simstudy/results/",dataset,"/plots/psm_fixed_npg_2025.png"),
     width     = 5,
     height    = 5,
     units     = "in",
     res       = 100,
     pointsize = 4)

# Display the 4x4 grid of plots
print(grid_plot)

dev.off()

# ((gg_psms[[1]] | gg_psms[[2]] | gg_psms[[3]] | gg_psms[[4]]) &
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))/
#   ( (gg_psms[[5]] | gg_psms[[6]] | gg_psms[[7]] | gg_psms[[8]]) &
#       theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))) &
#   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
# set.seed(3)
# N <- 300
# D = 2
# true_weights <- c(0.2, 0.2, 0.2, 0.3, 0.1)
# true_mus_1 <- c(-3, -3, 3, 3, -1)
# true_mus_2 <- c(-2.5, 3, -3, 3, 0)
# true_mus <- rbind(true_mus_1, true_mus_2)
# C <- length(true_mus_1)
#
# true_Sigmas <- array(
#   matrix(c(3, 1, 1, 3), 2, 2),
#   dim = c(2, 2, C)
# )
# true_Sigmas[,,5] <- matrix(c(0.25, 0, 0, 0.25), 2, 2)
#
# true_allocs <- sample(C, N, prob = true_weights, replace = TRUE)
#
# y <- matrix(NA, N, D)
# for(i in 1:N){
#
#   y[i,] <- rmvnorm(1, mean = true_mus[,true_allocs[i]], sigma = true_Sigmas[,,true_allocs[i]])
#
# }
#
# plot(y)

##############################
# M_a
##############################

# Fixed values only

# filtering out values depending on needs

pdf(paste0("simstudy/results/",dataset,"/plots/M_a_bar_fixed_2025.pdf"), width = 10, height = 10)

gg_M_a_bar_fixed <-
all_M_a.w %>%
  filter(!(gamma %in% c(0.5,99:101)) & !(zeta %in% c(98, 1, 3))) %>% arrange(gamma, zeta) %>%
  #filter(!(gamma %in% c(0.5,99:101)) & !(zeta %in% c(1, 3, 98))) %>%
  #filter(!(gamma %in% c(0.5,99:101))) %>%
  mutate(gamma = ifelse(gamma == 98, "G(3,2)", gamma)) %>%
  # mutate(gamma = ifelse(gamma == 100, "rp", gamma)) %>%
  # mutate(gamma = ifelse(gamma == 101, "rf", gamma)) %>%
  mutate(zeta = ifelse(zeta == 98, "G(3,2)", zeta)) %>%
  mutate(zeta = ifelse(zeta == 99, "G(0.1,1)", zeta)) %>%
  ggplot(aes(x = M_a, fill = factor(gamma)))+
  geom_bar(aes(y = after_stat(prop)),position = "dodge2") +
  facet_grid(gamma ~zeta, scales = "fixed") +
  labs(title = "",
       x = "",
       y = "",
       fill = expression(gamma)) +
  theme_minimal() +
  scale_fill_npg() +
  #theme(legend.position = "bottom") +
  theme(panel.spacing = unit(1, "cm"),  # Increase the distance between plots
    text=element_text(size=20), #change font size of all text
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        axis.text=element_text(size=15), #change font size of axis text
        axis.title=element_text(size=15), #change font size of axis titles
        plot.title=element_text(size=15), #change font size of plot title
        legend.text=element_text(size=15), #change font size of legend text
        legend.title=element_text(size=15)) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 10))
  #scale_x_continuous(breaks = seq(floor(min(all_M_a.w$M_a)), ceiling(max(all_M_a.w$M_a)), by = 1))

print(gg_M_a_bar_fixed)

dev.off()


# pdf(paste0("simstudy/results/",dataset,"/plots/M_a_bar_ratio.pdf"), width = 8, height = 4)
#
# gg_M_a_bar_ratio <-
# all_M_a.w %>%
#   filter(gamma %in% c(100) ) %>%
#   filter(!zeta %in% c(98)  ) %>%
#   mutate(gamma = ifelse(gamma == 98, "prior", gamma)) %>%
#   mutate(gamma = ifelse(gamma == 100, "", gamma)) %>%
#   mutate(gamma = ifelse(gamma == 101, "rf", gamma)) %>%
#   mutate(zeta = ifelse(zeta == 98, 100, zeta)) %>%
#   filter(gamma != "rf") %>%
#   ggplot(aes(x = M_a, fill = factor(zeta)))+
#   geom_bar(aes(y = after_stat(prop)),position = "dodge2") +
#   facet_grid(gamma ~zeta, scales = "fixed") +
#   labs(title = "",
#        x = "",
#        y = "",
#        fill = expression(zeta /gamma)) +
#   theme_minimal() +
#   scale_fill_npg() +
#   theme(legend.position = "bottom") +
#   theme(panel.spacing = unit(1, "cm"),
#     text=element_text(size=20), #change font size of all text
#         #strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         axis.text=element_text(size=15), #change font size of axis text
#         axis.title=element_text(size=15), #change font size of axis titles
#         plot.title=element_text(size=15), #change font size of plot title
#         legend.text=element_text(size=15), #change font size of legend text
#         legend.title=element_text(size=15))+
#   scale_x_continuous(breaks = scales::breaks_pretty(n = 10))
#
# print(gg_M_a_bar_ratio)
#
# dev.off()

#gg_M_a_bar_fixed / gg_M_a_bar_ratio


##########################
# PSM (Old color)
##########################


# pdf(paste0("simstudy/results/",dataset,"/plots/psm.pdf"), width = 15, height = 10)

# all_alphas <- c( 1)
# all_gammas <- c(0, 0.25, 0.5,  1, 98, 100, 101)
# all_zetas <- c(0.001,  1,  3, 5, 98)
#
# repulsive_grid_ordered <- expand.grid(all_alphas, all_zetas, all_gammas)
#
# layout(matrix(c(1:nrow(repulsive_grid)), length(unique(repulsive_grid$gamma)), length(unique(repulsive_grid$zeta)),byrow = TRUE))
# for(run in 1:nrow(repulsive_grid)){
#
#   filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
#   allocs <- readRDS(paste0("simstudy/results/",dataset,"/allocs/","allocs_",filename,".rds"))
#   #plotpsm(psm(allocs), main = paste0("PSM for ", gsub("_", " = ", filename)), method = "complete")
#   psm <- psm(allocs)
#   hc = hclust(as.dist(1 - psm), method = "complete", members = NULL)
#   psm_hc = psm
#   n = nrow(psm)
#   psm_hc[1:n, ] = psm_hc[hc$order, ]
#   psm_hc[, 1:n] = psm_hc[, hc$order]
#   image(1:n, 1:n, 1 - psm_hc, col = viridis(20))
#
#
# }

# Split them into fixed and prior values:
# Fixed:
repulsive_grid_fixed <- data.frame(repulsive_grid) %>%  mutate(row = row_number()) %>%
  filter(!(gamma %in% c(0.5,99:101)) & !(zeta %in% c(98:101)))


psm_fixed_plot_order <- repulsive_grid_fixed[order(repulsive_grid_fixed[, 2],repulsive_grid_fixed[, 3],repulsive_grid_fixed[, 1]), ] %>%
  pull(row)

#pdf(paste0("simstudy/results/",dataset,"/plots/psm_fixed.pdf"), width = 10, height = 10)

png( paste0("simstudy/results/",dataset,"/plots/psm_fixed.png"),
     width     = 5,
     height    = 5,
     units     = "in",
     res       = 500,
     pointsize = 4)

layout(matrix(c(1:nrow(repulsive_grid_fixed)), length(unique(repulsive_grid_fixed$gamma)), length(unique(repulsive_grid_fixed$zeta)),byrow = TRUE))
par(mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))
for(run in psm_fixed_plot_order){

  filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
  allocs <- readRDS(paste0("simstudy/results/",dataset,"/allocs/","allocs_",filename,".rds"))
  #plotpsm(psm(allocs), main = paste0("PSM for ", gsub("_", " = ", filename)), method = "complete")
  psm <- psm(allocs)
  hc = hclust(as.dist(1 - psm), method = "complete", members = NULL)
  psm_hc = psm
  n = nrow(psm)
  psm_hc[1:n, ] = psm_hc[hc$order, ]
  psm_hc[, 1:n] = psm_hc[, hc$order]
  image(1:n, 1:n, 1 - psm_hc, col = viridis(20), xlab = "", ylab = "")


}
dev.off()

# prior:
repulsive_grid_prior <- data.frame(repulsive_grid) %>%  mutate(row = row_number()) %>%
  filter((gamma %in% c(100))) %>%
  filter(!zeta%in% c(98)  )

psm_prior_plot_order <-  repulsive_grid_prior[order(repulsive_grid_prior[, 2],repulsive_grid_prior[, 3],repulsive_grid_prior[, 1]), ] %>%
  pull(row)

#switch the last 6 since ggplot switched rf and rp:
#psm_prior_plot_order2 <- c(psm_prior_plot_order[1:3],  psm_prior_plot_order[7:9], psm_prior_plot_order[4:6])

#pdf(paste0("simstudy/results/",dataset,"/plots/psm_prior.pdf"), width = 10, height = 10)

# png( paste0("simstudy/results/",dataset,"/plots/psm_ratio.png"),
#      width     = 8,
#      height    = 4,
#      units     = "in",
#      res       = 500,
#      pointsize = 4)

png( paste0("simstudy/results/",dataset,"/plots/psm_ratio_flat.png"),
     width     = 8,
     height    = 2,
     units     = "in",
     res       = 500,
     pointsize = 4)

layout(matrix(c(1:nrow(repulsive_grid_prior)), length(unique(repulsive_grid_prior$gamma)), length(unique(repulsive_grid_prior$zeta)),byrow = TRUE))
par(mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))
for(run in psm_prior_plot_order){

  filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
  allocs <- readRDS(paste0("simstudy/results/",dataset,"/allocs/","allocs_",filename,".rds"))
  #plotpsm(psm(allocs), main = paste0("PSM for ", gsub("_", " = ", filename)), method = "complete")
  psm <- psm(allocs)
  hc = hclust(as.dist(1 - psm), method = "complete", members = NULL)
  psm_hc = psm
  n = nrow(psm)
  psm_hc[1:n, ] = psm_hc[hc$order, ]
  psm_hc[, 1:n] = psm_hc[, hc$order]
  image(1:n, 1:n, 1 - psm_hc, col = viridis(20), xlab = "", ylab = "")


}

dev.off()


# Binder estimate:
#run = 1

pdf(paste0("simstudy/results/",dataset,"/plots/binder.pdf"), width = 12, height = 10)

for(start_index in seq(1, nrow(repulsive_grid), by = length(unique(all_gammas)))){
  do.call(grid.arrange, c(c(binder.ggs[start_index:(start_index + 3)]), ncol=2))
}

dev.off()


pdf(paste0("simstudy/results/",dataset,"/plots/vi.pdf"), width = 12, height = 10)

for(start_index in seq(1, nrow(repulsive_grid), by = length(unique(all_gammas)))){
  do.call(grid.arrange, c(c(vi.ggs[start_index:(start_index + 3)]), ncol=2))
}

dev.off()



################################################################################
# Gamma posteriors for the ratio plots
################################################################################

repulsive_grid_prior <- repulsive_grid %>% filter((gamma %in% c(100)) & (zeta %in% c(0.001, 0.1, 5, 99))) %>% arrange(gamma, zeta)



all_repara.w <- data.frame()
for(run in 1:nrow(repulsive_grid_prior)){

  filename <- paste0( paste0(colnames(repulsive_grid_prior[run,]), "_"), paste0(repulsive_grid_prior[run,]), collapse = "_" )

  M_a.w <- cbind(repulsive_grid_prior[run,],
                 M_a = readRDS(paste0("simstudy/results/",dataset,"/M_a/","M_a_",filename,".rds")))%>%
    mutate(iter = 1:n_save)


  repara.w <- cbind(repulsive_grid_prior[run,],
                    gamma_out = readRDS(paste0("simstudy/results/",dataset,"/gammas/gammas_",filename,".rds")),
                    zeta = readRDS(paste0("simstudy/results/",dataset,"/zetas/zetas_",filename,".rds")))%>%
                    mutate(iter = 1:n_save)


  all_repara.w <- rbind(all_repara.w, repara.w)


}

all_repara.l <- all_repara.w %>%
  pivot_longer(cols =-c("alpha", "gamma", "zeta", "iter"), values_to = "value", names_to = "repara")



# create_paper_plots <- function(df, geom, ncol, legend){
#
#
#   if(geom == "line"){
#
#     plot <-   df %>%
#       filter(repara == "gamma_out") %>%
#       mutate(repara = ifelse(repara == "gamma_out", "", repara)) %>%
#       filter(zeta != 3) %>%
#       mutate(zeta = ifelse(zeta == 99, 100, zeta)) %>%
#       ggplot(aes(x = iter, y = value, col = factor(zeta))) +
#       geom_line()+
#       facet_wrap(gamma ~ zeta, scales = "free", ncol = ncol) +
#       labs(title = "",
#            x = "",
#            y = "",
#            col = expression(gamma)) +
#       theme_minimal() +
#       scale_color_npg() +
#       if(legend == TRUE){
#         theme(legend.position = "bottom")
#       }else{
#         theme(legend.position = "none")
#       }  +
#       theme(panel.spacing = unit(1, "cm"),  # Increase the distance between plots
#             text=element_text(size=20), #change font size of all text
#             strip.text.x = element_blank(),
#             axis.text=element_text(size=15), #change font size of axis text
#             axis.title=element_text(size=15), #change font size of axis titles
#             plot.title=element_text(size=15), #change font size of plot title
#             legend.text=element_text(size=15), #change font size of legend text
#             legend.title=element_text(size=15))
#
#   }else if(geom == "density"){
#
#     plot <- df %>%
#       filter(repara == "gamma_out") %>%
#       filter(zeta != 3) %>%
#       mutate(zeta = ifelse(zeta == 99, 100, zeta)) %>%
#       ggplot(aes(x = value, fill = factor(zeta), col = factor(zeta))) +
#       geom_density(alpha = 0.4)+
#       facet_wrap(zeta ~ ., scales = "free", ncol = ncol) +
#       labs(title = "",
#            x = "",
#            y = "",
#            col = expression(gamma),
#            fill = expression(gamma)) +
#       theme_minimal() +
#       scale_fill_npg() +
#       scale_color_npg() +
#       if(legend == TRUE){
#         theme(legend.position = "bottom")
#       }else{
#         theme(legend.position = "none")
#       }  +
#       theme(panel.spacing = unit(1, "cm"),  # Increase the distance between plots
#             text=element_text(size=20), #change font size of all text
#             #strip.background = element_blank(),
#             strip.text.x = element_blank(),
#             axis.text=element_text(size=15), #change font size of axis text
#             axis.title=element_text(size=15), #change font size of axis titles
#             plot.title=element_text(size=15), #change font size of plot title
#             legend.text=element_text(size=15), #change font size of legend text
#             legend.title=element_text(size=15))
#
#
#   }else if(geom == "barplot"){
#
#     plot <-  df %>%
#       filter(zeta != 3) %>%
#       mutate(zeta = ifelse(zeta == 99, 100, zeta)) %>%
#       ggplot(aes(x = M_a, fill = factor(zeta)))+
#       geom_bar(aes(y = after_stat(prop)),position = "dodge2") +
#       facet_wrap(zeta ~., scales = "fixed", ncol = ncol) +
#       labs(title = "",
#            x = "",
#            y = "",
#            fill = expression(zeta /gamma)) +
#       theme_minimal() +
#       scale_fill_npg() +
#       # if(legend == TRUE){
#       #   theme(legend.position = "bottom")
#       # }else{
#       #   theme(legend.position = "none")
#       # }  +
#       theme(panel.spacing = unit(1, "cm"),
#             text=element_text(size=20), #change font size of all text
#             strip.background = element_blank(),
#             strip.text.x = element_blank(),
#             axis.text=element_text(size=15), #change font size of axis text
#             axis.title=element_text(size=15), #change font size of axis titles
#             plot.title=element_text(size=15), #change font size of plot title
#             legend.text=element_text(size=15), #change font size of legend text
#             legend.title=element_text(size=15))+
#       scale_x_continuous(breaks = scales::breaks_pretty(n = 10))
#
#   }
#
#   return(plot)
#
# }

# New function:

# to delete
# test <- all_repara.l%>%
#   filter(!(gamma %in% gamma_exclude) & !(zeta %in% zeta_exclude))
# table(test$gamma)
# table(test$zeta)

#create_paper_plots(all_repara.l, "line", 1,legend = FALSE, gamma_include = c(98, 101), zeta_include = c(0.001, 0.1, 5, 99))

#df = all_repara.l

repulsive_grid_prior

create_paper_plots <- function(df, geom, ncol, legend, gamma_include = c(98, 101), zeta_include = c(98)) {

  # Apply filtering based on the values to include
  if (!is.null(gamma_include)) {
    df <- df %>% filter(gamma %in% gamma_include)
  }

  if (!is.null(zeta_include)) {
    df <- df %>% filter(zeta %in% zeta_include)
  }


  df <- df %>%
    mutate(zeta = ifelse(zeta == 99, 100, zeta))


  if (geom == "line") {

    plot <- df %>%
      filter(repara == "gamma_out") %>%
      mutate(repara = ifelse(repara == "gamma_out", "", repara)) %>%
      ggplot(aes(x = iter, y = value, col = factor(zeta))) +
      geom_line() +
      facet_wrap(gamma ~ zeta, scales = "free", ncol = ncol) +
      labs(title = "", x = "", y = "", col = expression(gamma)) +
      theme_minimal() +
      scale_color_npg() +
      theme(panel.spacing = unit(1, "cm"),
            text = element_text(size = 20),
            strip.text.x = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            plot.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)) +
      if (legend) theme(legend.position = "bottom") else theme(legend.position = "none")

  } else if (geom == "density") {

    plot <- df %>%
      filter(repara == "gamma_out") %>%
      ggplot(aes(x = value, fill = factor(zeta), col = factor(zeta))) +
      geom_density(alpha = 0.4) +
      facet_wrap(zeta ~ ., scales = "free", ncol = ncol) +
      labs(title = "", x = "", y = "", col = expression(gamma), fill = expression(gamma)) +
      theme_minimal() +
      scale_fill_npg() +
      scale_color_npg() +
      theme(panel.spacing = unit(1, "cm"),
            text = element_text(size = 20),
            strip.text.x = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            plot.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)) +
      if (legend) theme(legend.position = "bottom") else theme(legend.position = "none")

  } else if (geom == "barplot") {

    plot <- df %>%
      ggplot(aes(x = M_a, fill = factor(zeta))) +
      geom_bar(aes(y = after_stat(prop)), position = "dodge2") +
      facet_wrap(zeta ~ ., scales = "fixed", ncol = ncol) +
      labs(title = "", x = "", y = "", fill = expression(zeta / gamma)) +
      theme_minimal() +
      scale_fill_npg() +
      theme(panel.spacing = unit(1, "cm"),
            text = element_text(size = 20),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            plot.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)) +
      scale_x_continuous(breaks = scales::breaks_pretty(n = 10))
  }

  return(plot)
}

path <- paste0("simstudy/results/",dataset,"/plots/")

# Plots and psm:

#pdf( paste0(path,"psm_ratio_ncol1_2.png"), width = 8, height = 16)

repulsive_grid_prior

gg_psms <- list()
for(run in 1:nrow(repulsive_grid_prior)){

  filename <- paste0( paste0(colnames(repulsive_grid_prior[run,]), "_"), paste0(repulsive_grid_prior[run,]), collapse = "_" )
  allocs <- readRDS(paste0("simstudy/results/",dataset,"/allocs/allocs_",filename,".rds"))

  #ind <- which(run == psm_prior_plot_order)

  # Compute the Posterior Similarity Matrix from the allocation matrix
  posterior_similarity <- matrix(0, nrow = ncol(allocs), ncol = ncol(allocs))

  for (i in 1:nrow(allocs)) {
    posterior_similarity <- posterior_similarity + (outer(allocs[i, ], allocs[i, ], "=="))
  }

  # Normalize the similarity matrix by the number of iterations
  posterior_similarity <- posterior_similarity / nrow(allocs)

  # Convert the matrix to a format suitable for ggplot2
  library(reshape2)
  similarity_df <- melt(posterior_similarity)
  colnames(similarity_df) <- c("Observation1", "Observation2", "Similarity")

  # Perform hierarchical clustering on the similarity matrix
  hc <- hclust(as.dist(1 - posterior_similarity))
  reordered_indices <- hc$order
  similarity_df$Observation1 <- factor(similarity_df$Observation1, levels = reordered_indices)
  similarity_df$Observation2 <- factor(similarity_df$Observation2, levels = reordered_indices)

  # Very Light Gray: #e0e0e0
  # Ivory: #fffff0
  # Lavender: #f4f3ff
  # Beige: #f5f5dc
  # Mint Cream: #f5fff5


  # good: "#f0f0f0"

  # Plot the clustered heatmap
  gg_psms[[run]] <- ggplot(similarity_df, aes(x = Observation1, y = Observation2, fill = Similarity)) +
    geom_tile() +
    scale_fill_gradient(low = "#f4f3ff", high = pal_npg()(4)[run]) +
    #scale_fill_gradient(low = "white", high = "darkred") +
    labs(title = "",
         x = "",
         y = "",
         fill = "") +
    theme_minimal() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

}

gg_psm_nps <- gg_psms[[1]] / gg_psms[[2]] / gg_psms[[3]] / gg_psms[[4]]

# Ncol 1 combined now with psm

# Example ggplot objects
p1 <- create_paper_plots(all_repara.l, "line", 1,legend = FALSE, gamma_include = c(100), zeta_include = c(0.001, 0.1, 5, 99))
p2 <-  create_paper_plots(all_repara.l, "density", 1,legend = FALSE, gamma_include = c(100), zeta_include = c(0.001, 0.1, 5, 99))
p3 <- create_paper_plots(all_M_a.w, "barplot", 1 ,legend = TRUE, gamma_include = c(100), zeta_include = c(0.001, 0.1, 5, 99)) +
  theme(
  legend.text = element_text(size = 20),  # Adjust legend text size
  legend.title = element_text(size = 20), # Adjust legend title size
  legend.key.size = unit(1.5, "cm")       # Adjust legend key size
)


#pdf(paste0(path,"M_bar_repara_line_dens_psm_legendside.pdf"), width = 32, height = 16)

png( paste0(path,"M_bar_repara_line_dens_psm_legendside_2025.png"),
     width     = 32,
     height    = 16,
     units     = "in",
     res       = 100,
     pointsize = 4)

# Combine plots with a shared legend

#(p1 + theme(legend.position = "none"))  | (p2 + theme(legend.position = "none"))| (gg_psm_nps + theme(legend.position = "none")) | (p3 + plot_layout(guides = "collect"))


p1 + p2 + gg_psm_nps +  p3 +
  plot_layout(guides = "collect", ncol = 4, widths = c(2, 1, 1, 2))

dev.off()



































































#################################################################################
# Old
#################################################################################





if(plots){
  dir.create(paste0("simstudy/results/",dataset,"/plots/"), showWarnings = FALSE, recursive = TRUE)
  #library(label.switching)
  library(mcclust.ext)
  library(reshape2)
  all_mus.w <- data.frame()
  all_weights.w <- data.frame()
  all_sigmas.w <- data.frame()
  all_M_a.w <- data.frame()
  all_post_dens.w <- data.frame()
  all_psms <- list()
  binder.ggs <- list()
  vi.ggs <- list()

  for(run in 1:nrow(repulsive_grid)){

    filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )

    # M_a

    M_a.w <- cbind(repulsive_grid[run,],
                   readRDS(paste0("simstudy/results/",dataset,"/M_a/","M_a_",filename,".rds")))%>%
      mutate(iter = 1:n_save)

    all_M_a.w <- rbind(all_M_a.w, M_a.w)


    # post_dens.w <- cbind(repulsive_grid[run,],y_line,
    #                      readRDS(paste0("simstudy/results/",dataset,"/post_dens/","post_dens_",filename,".rds")))
    #
    # all_post_dens.w <- rbind(all_post_dens.w, post_dens.w)

    zetas_samp <- rep(repulsive_grid$zeta[run], 2)
    gamma_samp <- repulsive_grid$gamma[run]
    alpha_prior <- repulsive_grid$alpha[run]
    filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )

    allocs <- readRDS(paste0("simstudy/results/",dataset,"/allocs/","allocs_",filename,".rds"))
    #allocs <- results$alloc_out
    s_binder <- c(salso(allocs, loss=salso::binder(a = 1)))


    binder.ggs[[run]] <-
      ggplot() +
      geom_point(data = df, aes(x = x, y = y, color = Cluster),  size = 2) +
      geom_encircle(data = df, aes(x = x, y = y, group = Cluster, fill = Cluster),
                    size = 0.5, expand = 0, alpha = 0.15, s_shape = 0.5,  spread=0.025) +
      labs(title = paste0("Binder \n zeta = ", zetas_samp[1], " and gamma = ", gamma_samp), x = "BMI", y = "CEBQ Person Parameter", legend = "test") +
      theme_minimal()

    s_vi <- c(salso(allocs, loss=salso::VI(a = 1)))

    df <- data.frame(x = y[,1], y = y[,2], Cluster = factor(s_vi))

    vi.ggs[[run]] <-
      ggplot() +
      geom_point(data = df, aes(x = x, y = y, color = Cluster),  size = 2) +
      geom_encircle(data = df, aes(x = x, y = y, group = Cluster, fill = Cluster),
                    size = 0.5, expand = 0, alpha = 0.15, s_shape=0.5,  spread=0.025) +
      labs(title = paste0("vi \n zeta = ", zetas_samp[1], " and gamma = ", gamma_samp), x = "BMI", y = "CEBQ Person Parameter", legend = "test") +
      theme_minimal()


  }


  # Binder estimate:
  #run = 1

  pdf(paste0("simstudy/results/",dataset,"/plots/binder.pdf"), width = 12, height = 10)

  for(start_index in seq(1, nrow(repulsive_grid), by = length(unique(all_gammas)))){
    do.call(grid.arrange, c(c(binder.ggs[start_index:(start_index + 3)]), ncol=2))
  }

  dev.off()


  pdf(paste0("simstudy/results/",dataset,"/plots/vi.pdf"), width = 12, height = 10)

  for(start_index in seq(1, nrow(repulsive_grid), by = length(unique(all_gammas)))){
    do.call(grid.arrange, c(c(vi.ggs[start_index:(start_index + 3)]), ncol=2))
  }

  dev.off()



  grid_rev <- data.frame(alpha = expand.grid(all_gammas, all_alphas)[,2], gamma = expand.grid(all_gammas, all_alphas)[,1])
  all_clust_vars.l <- data.frame()
  all_clust_ginis.l <- data.frame()
  all_clust_ents.l <- data.frame()

  for(run in 1:nrow(grid_rev)){
    print(paste0("run = ", run))
    filename <- paste0( paste0(colnames(grid_rev[run,]), "_"), paste0(grid_rev[run,]), collapse = "_" )

    # allocs
    allocs <- readRDS(paste0("simstudy/results/",dataset,"/allocs/","allocs_",filename,".rds"))

    clust_vars <- cbind(grid_rev[run,], value = variance_cluster_sizes(allocs))%>% mutate(iter = 1:n_save)
    all_clust_vars.l <- rbind(all_clust_vars.l, clust_vars)

    clust_ginis <- cbind(grid_rev[run,], value = gini_coefficient(allocs))%>% mutate(iter = 1:n_save)
    all_clust_ginis.l <- rbind(all_clust_ginis.l, clust_ginis)

    clust_ents <- cbind(grid_rev[run,], value = entropy(allocs))%>% mutate(iter = 1:n_save)
    all_clust_ents.l <- rbind(all_clust_ents.l, clust_ents)

  }

  #layout(matrix(c(1:nrow(repulsive_grid)), length(all_alphas), length(all_gammas), byrow = TRUE))

  grid_rev_small <- data.frame(alpha = expand.grid(all_gammas, all_alphas)[,2], gamma = expand.grid(all_gammas, all_alphas)[,1]) %>%
    filter(!gamma %in% c(0.5, 5))



  # pdf(paste0("simstudy/results/",dataset,"/plots/psm.pdf"), width = 15, height = 10)
  #
  layout(matrix(c(1:nrow(repulsive_grid)), length(unique(repulsive_grid$zeta)), length(unique(repulsive_grid$gamma)), byrow = TRUE))
  for(run in 1:nrow(repulsive_grid)){

    filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
    allocs <- readRDS(paste0("simstudy/results/",dataset,"/allocs/","allocs_",filename,".rds"))
    #plotpsm(psm(allocs), main = paste0("PSM for ", gsub("_", " = ", filename)), method = "complete")
    psm <- psm(allocs)
    hc = hclust(as.dist(1 - psm), method = "complete", members = NULL)
    psm_hc = psm
    n = nrow(psm)
    psm_hc[1:n, ] = psm_hc[hc$order, ]
    psm_hc[, 1:n] = psm_hc[, hc$order]
    image(1:n, 1:n, 1 - psm_hc, col = viridis(20))


  }
  #
  # dev.off()

  # New code for psm with ggplot:

  # pdf(paste0("simstudy/results/",dataset,"/plots/psm.pdf"), width = 15, height = 10)
  #

  plotpsm(psm(allocs), main = paste0("PSM for ", gsub("_", " = ", filename)), method = "complete")


  all_psm <- data.frame()
  gg_psm <- list()
  for(run in 1:nrow(repulsive_grid)){

    filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
    allocs <- readRDS(paste0("simstudy/results/",dataset,"/allocs/","allocs_",filename,".rds"))
    # Compute the Posterior Similarity Matrix from the allocation matrix
    psm <- matrix(0, nrow = ncol(allocs), ncol = ncol(allocs))

    for (i in 1:nrow(allocs)) {
      psm <- psm + (outer(allocs[i, ], allocs[i, ], "=="))
    }

    # Normalize the similarity matrix by the number of iterations
    psm <- psm / nrow(allocs)

    # Convert the matrix to a format suitable for ggplot2

    all_psm <- melt(posterior_similarity) %>% rbind(all_psm)

    # Perform hierarchical clustering on the similarity matrix
    #"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
    hc <- hclust(as.dist(1 - psm))
    reordered_indices <- hc$order
    similarity_df$Observation1 <- factor(similarity_df$Observation1, levels = reordered_indices)
    similarity_df$Observation2 <- factor(similarity_df$Observation2, levels = reordered_indices)

    hc = hclust(as.dist(1 - psm), method = method, members = NULL)
    psm_hc = psm
    n = nrow(psm)
    psm_hc[1:n, ] = psm_hc[hc$order, ]
    psm_hc[, 1:n] = psm_hc[, hc$order]
    image(1:n, 1:n, 1 - psm_hc, col = heat.colors(20), ...)

    psm_hc2 <- psm
    psm_hc2[1:n, ] = psm_hc2[reordered_indices,]
    psm_hc2[, 1:n] = psm_hc2[,reordered_indices]
    image(1:n, 1:n, 1 - psm_hc2, col = viridis(20))

    # Plot the clustered heatmap
    gg_psm[[run]] <- ggplot(similarity_df, aes(x = Observation1, y = Observation2, fill = Similarity)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "darkred") +
      labs(title = "",
           x = "",
           y = "",
           fill = "") +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank())

    ggplot(similarity_df, aes(x = Observation1, y = Observation2, fill = 1 - Similarity)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "darkred") +
      labs(title = "",
           x = "",
           y = "",
           fill = "") +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank())


  }

  colnames(all_psm) <- c("Observation1", "Observation2", "Similarity")


  library(gridExtra)
  grid.arrange(grobs = gg_psm, ncol = 4, nrow = 4)





  #



















  gini_var_ent_means <- rbind(
    all_clust_ents.l %>%
      dplyr::group_by(alpha, gamma) %>%
      dplyr::summarize(mean = mean ( value)) %>%
      mutate(variable = "Entropy")
    ,
    all_clust_vars.l %>%
      dplyr::group_by(alpha, gamma) %>%
      dplyr::summarize(mean = mean ( value)) %>%
      mutate(variable = "Variance")
    ,
    all_clust_ginis.l %>%
      dplyr::group_by(alpha, gamma) %>%
      dplyr::summarize(mean = mean ( value))%>%
      mutate(variable = "Gini Coefficient")
  ) %>%
    mutate(alpha = as.factor(alpha), gamma = as.factor(gamma))

  gini_var_ent_means %>%
    ggplot(aes(x = alpha, y = mean, fill = gamma)) +
    geom_bar(stat = "identity", position = position_dodge2()) +
    facet_wrap(variable~., scales = "free", ncol = 1) +
    theme_minimal()

  pdf(paste0("simstudy/results/",dataset,"/plots/gini_var.pdf"), width = 15, height = 10)

  gini_var_ent_means %>%
    filter(variable!="Entropy", !gamma %in% c(0.5, 5)) %>%
    ggplot(aes(x = alpha, y = mean, fill = gamma)) +
    geom_bar(stat = "identity", position = position_dodge2()) +
    facet_wrap(variable~., scales = "free", ncol = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
          strip.text.x = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title  = element_text(size = 20))

  dev.off()



  # M_a
  #M_a.l <- all_M_a.w %>% pivot_longer(cols = -c("alpha", "gamma", "zeta","iter"), names_to = "component", values_to = "value")

  colnames(all_M_a.w)[4] <- "value"

  pdf(paste0("simstudy/results/",dataset,"/plots/M_bar.pdf"), width = 15, height = 10)

  all_M_a.w %>%
    dplyr::mutate(across(any_of(c(colnames(repulsive_grid), "component")), as.factor)) %>%
    filter(!gamma %in% c(0.5, 5)) %>%
    ggplot(aes(x = value, fill = gamma)) +
    geom_bar(aes(y = ..prop..), position = position_dodge2(width = 0.9, preserve = "single")) +
    facet_wrap( zeta ~ . , scales = "free") +
    theme_bw() +
    ggtitle(paste0("")) +
    xlab("")+
    ylab("") +
    theme(text=element_text(size=20), #change font size of all text
          axis.text=element_text(size=20), #change font size of axis text
          axis.title=element_text(size=20), #change font size of axis titles
          plot.title=element_text(size=20), #change font size of plot title
          legend.text=element_text(size=25), #change font size of legend text
          legend.title=element_text(size=25)) #change font size of legend title

  dev.off()

  M_a.l %>%
    dplyr::mutate(across(any_of(c(colnames(repulsive_grid), "component")), as.factor)) %>%
    ggplot(aes(x = iter, y = value, col = gamma)) +
    geom_line() +
    facet_wrap( alpha ~ gamma , scales = "fixed", ncol = length(all_gammas))

  # Posterior density

  colnames(all_post_dens.w)[4] <- "post_dens"

  pdf(paste0("simstudy/results/",dataset,"/plots/post_dens.pdf"), width = 15, height = 10)

  all_post_dens.w %>%
    filter(!gamma %in% c(0.5, 5)) %>%
    dplyr::mutate(across(any_of(c(colnames(repulsive_grid), "component")), as.factor)) %>%
    ggplot() +
    geom_line(aes(x = y_line, y = post_dens, col = gamma), size = 1) +
    #scale_size_identity() +
    #scale_linetype_identity() +
    geom_histogram(data = data.frame(value = y), aes(x = value, y = ..density..), alpha = .5, bins = 40) +
    facet_wrap(alpha~gamma, ncol =3) +
    theme_bw() +
    ggtitle(paste0("")) +
    xlab("")+
    ylab("") +
    theme( #strip.background = element_blank(),
      #strip.text.x = element_blank(),
      text=element_text(size=20), #change font size of all text
      axis.text=element_text(size=20), #change font size of axis text
      axis.title=element_text(size=20), #change font size of axis titles
      plot.title=element_text(size=20), #change font size of plot title
      legend.text=element_text(size=25), #change font size of legend text
      legend.title=element_text(size=25)) #change font size of legend title

  dev.off()




}


