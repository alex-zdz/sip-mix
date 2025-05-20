
rep_mix_main <- function(y, updates, hyperparameters, starting_values,
                           n_save, n_burn, n_thin){
  
  #hyperparameters = NULL
  #starting_values = NULL
  #updates = NULL
  
  if(TRUE){
    #starting_values = NULL
    #set.seed(1)
    niter <- n_burn + n_thin * n_save
    N <- ifelse(is.null(dim(y)), length(y), nrow(y))
    D = ifelse(is.null(dim(y)), 1, ncol(y))
    
    if(D == 1){
      from_to <- range(y)
      y_line <- seq(from_to[1], from_to[2], length.out = 5e1)
    }else if(D == 2){
      from_to <- apply(y, 2, range)
      y_line1 <- seq(from_to[1,1], from_to[2,1], length.out = 2e1)
      y_line2 <- seq(from_to[1,2], from_to[2,2], length.out = 2e1)
      y_grid <- as.matrix(expand.grid(y_line1, y_line2))
    }else{
      
      print("Posterior evaluation not finished for D>2")
    }
    
    alloc_out <- matrix(NA, n_save, N)
    w_rank_out <- matrix(NA, n_save, 4)
    
    y_pred_out <- matrix(NA, n_save, D) # change
    M_a_out <- rep(NA, n_save)
    M_na_out <- rep(NA, n_save)
    # count number of births and deaths
    n_bd_out  <- matrix(NA, niter, 3)
    n_cs_out <- matrix(0, n_save, N)
    # repulsion parameters
    gamma_out <- rep(NA, n_save)
    zetas_out <- matrix(NA, n_save, D) # change
    # Birth & Deat
    bd_out <- rep(NA, n_save)
    log_prior_mu_out <- rep(NA, n_save)
    log_lik_mu_out <- rep(NA, n_save)
    log_accept_death_out <- rep(0, n_save)
    log_accept_birth_out <- rep(0, n_save)
    log_accept_birth_M_out <- rep(0, n_save)
    log_accept_birth_weights_out <- rep(0, n_save)
    log_accept_birth_locs_out <- rep(0, n_save)
    log_accept_birth_bd_props_out <- rep(0, n_save)
    # lik and loglik
    lik_out <- rep(0, n_save)
    loglik_out <- rep(0, n_save)
    # pearson distance
    pearson_cs_dist_out <- rep(NA, n_save)
    # posterior density estimate
    if(D == 1){
      post_dens_out <- matrix(NA, n_save, length(y_line))
    }else if (D == 2){
      post_dens_out <- matrix(NA, n_save, nrow(y_grid)) # change
    }else{
      post_dens_out <- matrix(NA, n_save, 1)
    }
    
    # locations log pwd 
    loc_log_pwd_out <- rep(0, n_save)
    # GE
    mu_norm_consts_out = matrix(NA, n_save, D)
    mu_kernels_out = matrix(NA, n_save, D)
    mu_pwd_terms_out = matrix(NA, n_save, D)
    
    
    if (is.null(hyperparameters)) {
      
      hyperparameters <- list(
        mu_prior_mean = mean(y),
        mu_prior_var_scale = 1,
        mu_prior_var = 100,
        mu_prop_means = if(D == 1) mean(y) else colMeans(y),
        mu_prop_Sigma = if(D == 1) 1 else diag(rep(10, D)),
        Sigma_prior_mat = diag(D),
        Sigma_prior_df = D,
        sigma2_prior_shape = 3,
        sigma2_prior_rate = 1,
        alpha_prior = 0.3,
        beta_prior = 0.3,
        M_prior = 2,
        p_b = 0.5,
        rep_prior_shape = 3,
        rep_prior_rate = 2,
        g_to_z_ratio = 1,
        g_time_z_ratio = 10,
        repulsive_locations = TRUE,
        mu_independent = TRUE,
        weights_prior = "selberg"
      )
      
    }else{
      
      hyperparameter_names <- c(
        "mu_prior_mean",
        "mu_prior_var_scale",
        "mu_prior_var",
        "mu_prop_means",
        "mu_prop_Sigma",
        "Sigma_prior_mat",
        "Sigma_prior_df",
        "sigma2_prior_shape",
        "sigma2_prior_rate",
        "alpha_prior",
        "beta_prior",
        "M_prior",
        "p_b",
        "rep_prior_shape",
        "rep_prior_rate",
        "g_to_z_ratio",
        "g_time_z_ratio",
        "repulsive_locations",
        "mu_independent",
        "weights_prior",
        "scaling_factor"
      )
      
      all_hyperparameters_supplied <- names(hyperparameters) %in% hyperparameter_names
      if(all(all_hyperparameters_supplied)){
        # check values
        print("All hyperparameters supplied")
      }else{
        print(paste0("No value supplied for ", hyperparameter_names[!all_hyperparameters_supplied]))
      }
      
    }
    
    # Extract the variables from the hyperparameters list and assign them to individual variables
    list2env(hyperparameters, envir = .GlobalEnv)
    
    # Starting values
    
    if (is.null(starting_values)) {
      
      M_a_samp = max(2, rpois(1, M_prior))
      M_na_samp = max(2, rpois(1, M_prior))
      M_samp <- M_a_samp + M_na_samp
      allocs_samp = rep(c(1:M_a_samp), length.out = N)
      # change
      
      if(D == 1){
        mus_na_samp = rnorm(M_a_samp, mu_prop_means, mu_prop_Sigma)
        mus_a_samp = rnorm(M_na_samp, mu_prop_means, mu_prop_Sigma)
        mus_samp <-  c(mus_a_samp, mus_na_samp)
      }else{
        mus_na_samp = t(rmvnorm(M_na_samp, rep(mu_prop_means, D), diag(rep(mu_prop_var, D))))
        mus_a_samp = t(rmvnorm(M_a_samp, rep(mu_prop_means, D), diag(rep(mu_prop_var, D))))
        mus_samp <- cbind(mus_a_samp, mus_na_samp)
      }
      
      if(D == 1){
        Sigmas_a_samp = rep(1, M_a_samp)
        Sigmas_na_samp = rep(1, M_na_samp)
        Sigmas_samp = c(Sigmas_a_samp, Sigmas_na_samp)
      }else{
        Sigmas_a_samp = array(diag(D), dim = c(D, D, M_a_samp))
        Sigmas_na_samp = array(diag(D), dim = c(D, D, M_na_samp))
        Sigmas_samp <- array(c(Sigmas_a_samp, Sigmas_na_samp), dim = c(D, D, M_samp))
      }
      
      weights_samp <- c(brms::rdirichlet(1, rep(alpha_prior, M_samp)))
      gamma_samp = 1
      zetas_samp = rep(0.1, D)
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
      
    }else{
      
      starting_values_names <- c(
        "M_a_samp", "M_na_samp", "M_samp",
        "allocs_samp",
        "mus_a_samp", "mus_na_samp","mus_samp",
        "Sigmas_a_samp", "Sigmas_na_samp",
        "weights_samp",
        "gamma_samp", "zetas_samp",
        "tuning_sd_mu", "tuning_sd_gamma", "tuning_sd_zeta"
      )
      
      all_starting_values_supplied <- names(starting_values) %in% starting_values_names
      if(all(all_starting_values_supplied)){
        # check values
        print("All starting values supplied")
      }else{
        print(paste0("No value supplied for ", starting_values_names[!all_starting_values_supplied]))
      }
      
    }
    
    # Extract the variables from the priors list and assign them to individual variables
    list2env(starting_values, envir = .GlobalEnv)
    
    n_curr <- c(table(allocs_samp), rep(0,M_na_samp))
    
    if (is.null(updates)) {
      
      updates <- list(
        update_allocs = TRUE, update_gamma = TRUE, update_zetas = TRUE, 
        g_to_z = FALSE, g_time_z = FALSE,
        update_var = TRUE, update_mu = TRUE, update_weight = TRUE,
        update_M = TRUE, update_post_dens = TRUE
      )
      
    }else{
      
      updates_names <- c("update_allocs", "update_gamma", "update_zetas", "g_to_z", "g_time_z",
                         "update_var", "update_mu", "update_weight", "update_M", "update_post_dens")
      
      all_updates_supplied <- names(updates) %in% updates_names
      if(all(all_updates_supplied)){
        # check values
      }else{
        print(paste0("No value supplied for ", updates_names[!all_updates_supplied]))
      }
      
    }
    
    # Extract the variables from the priors list and assign them to individual variables
    list2env(updates, envir = .GlobalEnv)
    
    # Define mu and sigma out:
    mu_out = Sigma_out = weights_out <- if(D == 1){
      matrix(NA, n_save, M_samp)
    }else{
      vector("list", n_save) 
    }
    
    bd <- NA
    log_prior_mu <- NA
    log_lik_mu <- NA
    log_accept_death <- NA
    log_accept_birth <- NA
    log_accept_birth_M <- NA
    log_accept_birth_weights <- NA
    log_accept_birth_locs <- NA
    log_accept_birth_bd_props <- NA
    
    if(g_to_z | g_time_z){
      #update_gamma = TRUE
      update_zetas = FALSE
    }
    
    if(gamma_samp == 0 & update_gamma == TRUE){
      stop("Starting value of gamma cannot be 0 unless fixed.")
    }
    
    # if(update_allocs == FALSE & update_var == TRUE){
    #   sigma2_post_rate_out <- matrix(NA, n_save, M_a_samp)
    #   sigma2_post_shape_out <- matrix(NA, n_save, M_a_samp)
    # }
    
    #Set M_na to 0 if we do not wish to update the number of components.
    if(update_M == FALSE ){
      # M_na_samp = 0
      # M_samp <- M_a_samp
      # mus_na_samp <- numeric(0)
      # mus_samp <- mus_a_samp
      # Sigmas_na_samp = numeric(0)
      # Sigmas_samp = Sigmas_a_samp
      birth_death <- NA
      actual_birth <- NA
      actual_death <- NA
    }
    
  }
  
  ################################################################################
  
  for(iter in 1:niter){
    # for(iter in 1:100){
    
    if(D > 1 &!is.matrix(mus_a_samp)) stop("wrong dim")
    
    if(update_allocs){
      
      # 1.a Sample indicators from posterior proportional to current weight times density:
      if(D > 1){
        
        if(M_na_samp == 0){

          allocs_list <- sample_allocs_mv(N, M_samp, M_a_samp, M_na_samp,
                                          allocs_samp, weights_samp,
                                          y, mus_a_samp,
                                          Sigmas_a_samp, mus_na_samp = t(matrix(NA, 1, D)),
                                          Sigmas_na_samp = array(NA, dim = c(D, D, 1)), n_curr)

        }else{
          allocs_list <- sample_allocs_mv(N, M_samp, M_a_samp, M_na_samp,
                                          allocs_samp, weights_samp,
                                          y, mus_a_samp,
                                          Sigmas_a_samp, mus_na_samp,
                                          Sigmas_na_samp, n_curr)
        }

        allocs_samp <- allocs_list$allocs_samp
        n_curr <- allocs_list$n_curr
        
        # alloc_rates <- rep(0, M_samp)
        # for(i in 1:N){
        #   old_alloc <- allocs_samp[i]
        #   
        #   # Calculate allocation rates for allocated clusters
        #   for(c in 1:M_a_samp){
        #     # dmvnorm for multivariate normal distribution
        #     alloc_rates[c] <- weights_samp[c] * dmvnorm(y[i, ], mus_a_samp[, c], sigma = Sigmas_a_samp[,, c])
        #   }
        #   
        #   # Calculate allocation rates for non-allocated clusters
        #   if(M_na_samp != 0){
        #     for(c in 1:M_na_samp){
        #       alloc_rates[M_a_samp + c] <- weights_samp[M_a_samp + c] * dmvnorm(y[i, ], mus_na_samp[, c], sigma = Sigmas_na_samp[,, c])
        #     }
        #   }
        #   
        #   # Normalize the allocation rates to get allocation probabilities
        #   alloc_probs <- alloc_rates / sum(alloc_rates)
        #   
        #   # Sample a new allocation based on the probabilities
        #   allocs_samp[i] <- sample(M_samp, 1, prob = alloc_probs)
        #   
        #   # Update the counts for allocations
        #   n_curr[old_alloc] <- n_curr[old_alloc] - 1
        #   n_curr[allocs_samp[i]] <- n_curr[allocs_samp[i]] + 1
        # }
        
      }else{
        
        
        # allocs_list <- sample_allocs_uv(N, M_samp, M_a_samp, M_na_samp,
        #                                 allocs_samp, weights_samp,
        #                                 y, mus_a_samp,
        #                                 Sigmas_a_samp, mus_na_samp,
        #                                 Sigmas_na_samp, n_curr)
        
        
        alloc_rates <- rep(0, M_samp)
        for(i in 1:N){
          old_alloc <- allocs_samp[i]
          for(c in 1:M_a_samp){
            alloc_rates[c] <- weights_samp[c] * dnorm(y[i], mus_a_samp[c], sd = sqrt(Sigmas_a_samp[c]))
          }
          if(M_na_samp != 0){
            for(c in 1:M_na_samp){
              alloc_rates[M_a_samp + c] <- weights_samp[M_a_samp + c] * dnorm(y[i], mus_na_samp[c], sd = sqrt(Sigmas_na_samp[c]))
            }
          }
          if(sum(alloc_rates) == 0){
            alloc_probs <- rep(1 / M_samp, M_samp)
          }else{
            alloc_probs <- alloc_rates / sum(alloc_rates)
          }
          
          allocs_samp[i] <- sample(M_samp, 1, prob = alloc_probs)
          n_curr[old_alloc] <- n_curr[old_alloc] - 1
          n_curr[allocs_samp[i] ] <- n_curr[allocs_samp[i] ] + 1
        }
        
      }
      
      
     
      
      # Reorder elements
      new_alloc <- sort(unique(allocs_samp[(allocs_samp > M_a_samp)])) - M_a_samp
      empty_comp <- rev(which(n_curr[1:M_a_samp] == 0))
      if(length(empty_comp) != 0 ){
        
        #for(missing in rev(which(n_curr == 0))){
        for(missing in empty_comp){
          allocs_samp[allocs_samp > missing ] <- allocs_samp[allocs_samp > missing] - 1
          M_a_samp <- M_a_samp - 1 
          M_na_samp <- M_na_samp + 1
          n_curr <- n_curr[ - missing]
          n_curr <- c(n_curr, 0)
          weights_samp <- c(weights_samp[-missing], weights_samp[missing])
        }
        if(length(n_curr) != M_a_samp + M_na_samp) stop("numbers dont add up")
        # remove emptied components from allocated vectors and add them to nonallocated ones
        
        if( D == 1){
          mus_na_samp <- c(mus_na_samp, mus_a_samp[empty_comp])
          mus_a_samp <- mus_a_samp[-empty_comp]
          Sigmas_na_samp <- c(Sigmas_na_samp, Sigmas_a_samp[empty_comp])
          Sigmas_a_samp <- Sigmas_a_samp[-empty_comp]
        }else{
          mus_na_samp <- cbind(mus_na_samp, mus_a_samp[,empty_comp])
          mus_a_samp <- mus_a_samp[,-empty_comp, drop = FALSE]
          Sigmas_na_samp <- array(c(Sigmas_na_samp, Sigmas_a_samp[,,empty_comp]), dim=c(D, D, M_na_samp) )
          Sigmas_a_samp <- Sigmas_a_samp[,,-empty_comp, drop = FALSE]
        }
      }
      if(length(new_alloc) != 0){ #Adding na-components to allocated ones
        
        for(new in sort(new_alloc + M_a_samp, decreasing = FALSE)){
          allocs_samp[allocs_samp == new] <- M_a_samp + 1
          weights_samp <- append(weights_samp[-new], weights_samp[new], after = M_a_samp)
          M_a_samp <- M_a_samp + 1 
          M_na_samp <- M_na_samp - 1
        }
        n_curr <- c(table(allocs_samp), rep(0, M_na_samp))
        # Add newly filled na-components
        # 1.b Update allocated vectors in case a non-allocated component has been filled
        # 1.c Remove the newly filled components from the na-vectors
        
        if( D == 1){
          mus_a_samp <- c(mus_a_samp, mus_na_samp[new_alloc])
          Sigmas_a_samp <- c(Sigmas_a_samp, Sigmas_na_samp[new_alloc])
          mus_na_samp <- mus_na_samp[- new_alloc]
          Sigmas_na_samp <- Sigmas_na_samp[- new_alloc]
          mus_samp <- c(mus_a_samp, mus_na_samp)
        }else{
          mus_a_samp <- cbind(mus_a_samp, mus_na_samp[,new_alloc])
          mus_na_samp <- mus_na_samp[,-new_alloc, drop = FALSE]
          Sigmas_a_samp <- array(c(Sigmas_a_samp, Sigmas_na_samp[,,new_alloc]), dim=c(D, D, M_a_samp) )
          Sigmas_na_samp <- Sigmas_na_samp[,,-new_alloc, drop = FALSE]
          mus_samp <- cbind(mus_a_samp, mus_na_samp)
        }
        
      }
      
      M_samp <- M_a_samp + M_na_samp
      
    }
    
    if( D == 1){
      if(length(mus_samp) != M_samp)stop("error after allocs")
    }else{
      if(ncol(mus_samp) != M_samp)stop("error after allocs")
    }
    
    # 2. Update component specific variance sigma2
    if(update_var){
      if(D > 1){
        Sigma_post_mat <- array(NA, dim = c(D, D, M_a_samp))
        for(c in 1:M_a_samp){
          n_c <-  sum(allocs_samp == c)
          y_c <- y[allocs_samp == c,]
          Sigma_post_mat[,,c] = if(n_c > 1){
            t((y_c - matrix(mus_samp[,c], n_c, D, byrow = TRUE)))%*%(y_c - matrix(mus_samp[,c], n_c, D, byrow = TRUE)) + Sigma_prior_mat
          }else{
            y_c %*% t(y_c) + Sigma_prior_mat
          }
          Sigma_post_df = n_c + Sigma_prior_df
          Sigmas_a_samp[,,c] <- riwish(v = Sigma_post_df, S = Sigma_post_mat[,,c])
        }
        # 2.b Update non-allcoated component specific parameters by sampling from the prior 
        if(M_na_samp > 0){
          for(m in 1:M_na_samp){
            Sigmas_na_samp[,,m] <- riwish(v = Sigma_prior_df, S = Sigma_prior_mat)
          }
        }
        # Save the following for later sampling mu:
        log_lik <- c( rep(NA, M_a_samp), rep(0 , M_na_samp))
        for(c in 1:M_a_samp){
          n_c <-  sum(allocs_samp == c)
          log_lik[c] <- sum(log(dmvnorm_arma_mv(y[allocs_samp == c,, drop = FALSE], mean = mus_samp[,c], sigma = Sigmas_a_samp[,,c])))
        }    
        Sigmas_samp <- array(c(Sigmas_a_samp, Sigmas_na_samp), dim = c(D, D, M_samp))
      }else{
        # 1D:
        
        if(mu_independent){
          
          for(c in 1:M_a_samp){
            sum_y_mu_2 <- ifelse(sum(allocs_samp == c) != 0, sum((y[allocs_samp == c] - mus_samp[c])^2), 0)
            mean_y <- ifelse(sum(allocs_samp == c) != 0, mean(y[allocs_samp == c]), 0)
            n_c <- sum(allocs_samp == c)
            sigma2_post_rate <- sum_y_mu_2 / 2 + sigma2_prior_rate
            Sigmas_a_samp[c] <- 1 / rgamma(1, shape = n_c / 2 + sigma2_prior_shape, rate = sigma2_post_rate)
          }
          
          if(M_na_samp != 0){
            Sigmas_na_samp <- 1 / rgamma(M_na_samp, shape = sigma2_prior_shape, rate = sigma2_prior_rate)
          }
          
        }else{
          
          for(c in 1:M_a_samp){
            sum_y_mu_2 <- ifelse(sum(allocs_samp == c) != 0, sum((y[allocs_samp == c] - mus_samp[c])^2), 0)
            mean_y <- ifelse(sum(allocs_samp == c) != 0, mean(y[allocs_samp == c]), 0)
            n_c <- sum(allocs_samp == c)
            sigma2_post_rate <- 1 / 2 * (sum_y_mu_2 + ((mus_samp[c] - mu_prior_mean)^2 * mu_prior_var_scale)) + sigma2_prior_rate
            Sigmas_a_samp[c] <- 1 / rgamma(1, shape = (n_c  + 1) / 2 + sigma2_prior_shape, rate = sigma2_post_rate)
          }
          
          if(M_na_samp != 0){
            Sigmas_na_samp <- 1 / rgamma(M_na_samp, shape = sigma2_prior_shape, rate = sigma2_prior_rate)
          }
          
        }
        
        Sigmas_samp <- c(Sigmas_a_samp, Sigmas_na_samp)
        
      } 
      
    }
    
    mu_norm_consts <- rep(NA, D)
    mu_kernels <- rep(NA, D)
    mu_pwd_terms <- rep(NA, D)
    
    if(update_mu){
      if(repulsive_locations){
        if(D == 1){
          
          sum_y_mu_2s <- rep(NA, M_a_samp)
          for(c in 1:M_a_samp){
            sum_y_mu_2s[c] <- ifelse(sum(allocs_samp == c) != 0, sum((y[allocs_samp == c] - mus_a_samp[c])^2), 0)
          }
          
          # Likelihood part with 0s for non-allocated components
          log_ge_lik <- c((sum_y_mu_2s / (2 * Sigmas_a_samp)), rep(0, M_na_samp))
          # pairwise differences part of the current locations
          all_log_pwdiff <- apply(as.matrix(rdist(c(mus_samp))), 1, function(x) sum(zetas_samp * log(abs(x))[!is.infinite(log(abs(x)))]))
          log_demon <- - zetas_samp / 2 * mus_samp^2 + zetas_samp * all_log_pwdiff - log_ge_lik
          for(p in 1:M_samp){ 
            # We use a symmetric normal proposal which cancels in the acceptance rate
            x_tmp <- rnorm(1, mean = mus_samp[p] , sd = tuning_sd_mu)
            x_prop <- ifelse(x_tmp == 0, exp(-50), x_tmp)
            y_mu_prop_2 <- ifelse(sum(allocs_samp == p) != 0, sum((y[allocs_samp == p] - x_prop)^2), 0)
            log_numer <-  - zetas_samp / 2 * x_prop^2 + sum(zetas_samp * log(abs(cdist(x_prop, mus_samp[-p])))) -
              # likelihood part
              y_mu_prop_2 / (2 * Sigmas_samp[p])
            log_acrate <- min(log_numer - log_demon[p], 1)
            if(log(runif(1)) < log_acrate){
              mus_samp[p] <-  x_prop
              log_demon[p] <- log_numer
            }
          }
          # split locations 
          mus_a_samp <- mus_samp[1:M_a_samp]
          if(M_na_samp != 0){
            mus_na_samp <- mus_samp[(M_a_samp + 1):M_samp]
          }
        }else{
          for(d in 1:D){
            
            # infinite part only fix for the case if starting values are the same
            all_log_pwdiff <- apply(as.matrix(rdist(c(mus_samp[d,]))), 1, function(x) sum(zetas_samp[d] * log(abs(x))[!is.infinite(log(abs(x)))]))
            log_demon <- - zetas_samp[d] / 2 * mus_samp[d, ]^2 + zetas_samp[d] * all_log_pwdiff + log_lik
            
            for(c in 1:M_samp){ 
              # We use a symmetric normal proposal which cancels in the acceptance rate
              x_tmp <- rnorm(1, mean = mus_samp[d, c] , sd = tuning_sd_mu)
              x_prop <- ifelse(x_tmp == 0, exp(-50), x_tmp)
              mu_c_new <- mus_samp[,c]
              mu_c_new[d] <- x_prop
              
              if(c <= M_a_samp){
                log_lik_c_new <- sum(log(dmvnorm_arma_mv(y[allocs_samp == c, , drop = FALSE], mean = mu_c_new, sigma = Sigmas_samp[,,c])))
              }else{
                log_lik_c_new <- 0
              }
              
              log_numer <-  - zetas_samp[d] / 2 * x_prop^2 + sum(zetas_samp[d] * log(abs(cdist(x_prop, mus_samp[d, -c])))) + log_lik_c_new
              log_acrate <- min(log_numer - log_demon[c], 1)
              if(log(runif(1)) < log_acrate){
                mus_samp[d,c] <-  x_prop
                log_demon[c] <- log_numer
              }
            }
            # split locations 
            mus_a_samp <- mus_samp[,1:M_a_samp,  drop = FALSE]
            if(M_na_samp != 0){
              mus_na_samp <- mus_samp[,(M_a_samp + 1):M_samp,  drop = FALSE]
            }
            
            # Save the individual parts of the GE distribution
            mu_norm_consts[d] <- log_G(M_a_samp, zetas_samp[d])
            mu_kernels[d] <- - sum(zetas_samp[d] / 2 * mus_a_samp[d, ]^2)
            mu_pwd_terms[d] <- zetas_samp[d] * sum(log(c(rdist(mus_a_samp[d,]))))
            
          }
        }
        
      }else{
        if(D == 1){
          # So far non-repulsive only in 1D
          if(mu_independent){
            for(c in 1:M_samp){
              sum_y_mu_2 <- ifelse(sum(allocs_samp == c) != 0, sum((y[allocs_samp == c] - mus_samp[c])^2), 0)
              mean_y <- ifelse(sum(allocs_samp == c) != 0, mean(y[allocs_samp == c]), 0)
              n_c <- sum(allocs_samp == c)
              # 2.a first the component means
              #post_var <-  1 / ((n_c / Sigmas_samp[c]) + mu_prior_var_scale / Sigmas_samp[c])   
              post_var <-  Sigmas_samp[c] / (n_c + mu_prior_var_scale)   
              post_mean <- post_var * ((mean_y * n_c + mu_prior_mean * mu_prior_var_scale) / Sigmas_samp[c])
              mus_samp[c] <- rnorm(1, mean = post_mean, sd = sqrt(post_var))
            }
          }else{
            for(c in 1:M_samp){
              sum_y_mu_2 <- ifelse(sum(allocs_samp == c) != 0, sum((y[allocs_samp == c] - mus_samp[c])^2), 0)
              mean_y <- ifelse(sum(allocs_samp == c) != 0, mean(y[allocs_samp == c]), 0)
              n_c <- sum(allocs_samp == c)
              # 2.a first the component means
              #post_var <-  1 / ((n_c / Sigmas_samp[c]) + mu_prior_var_scale / Sigmas_samp[c])   
              post_var <-  1 / (n_c / Sigmas_samp[c] + 1 / mu_prior_var)   
              post_mean <- post_var * ((mean_y * n_c) / Sigmas_samp[c] + mu_prior_mean / mu_prior_var) 
              mus_samp[c] <- rnorm(1, mean = post_mean, sd = sqrt(post_var))
            }
          }
          
          # split locations 
          mus_a_samp <- mus_samp[1:M_a_samp]
          if(M_na_samp != 0){
            mus_na_samp <- mus_samp[(M_a_samp + 1):M_samp]
          }
          
        }else{
      
          # Loop through each cluster
          for (c in 1:M_samp) {
            # Get the number of points allocated to cluster c
            n_c <- sum(allocs_samp == c)
            
            # Extract the data points belonging to cluster c
            y_c <- y[allocs_samp == c, ,  drop = FALSE]
            
            # Calculate the sample mean for the current cluster
            if (n_c > 0) {
              y_bar_c <- colMeans(y_c)  # Mean of the allocated points (Dx1 vector)
            } else {
              y_bar_c <- rep(0, ncol(y))  # If no data is assigned to the cluster, mean is zero (Dx1 vector)
            }
            
            # Posterior covariance for cluster c (DxD matrix)
            posterior_cov_c <- solve(n_c * solve(Sigmas_samp[,,c]) + solve(mu_prop_Sigma))
            
            # Posterior mean for cluster c (Dx1 vector)
            posterior_mean_c <- posterior_cov_c %*% (n_c * solve(Sigmas_samp[,,c]) %*% y_bar_c + solve(mu_prop_Sigma) %*% mu_prop_means)
            
            mus_samp[,c] <- rmvnorm(1, posterior_mean_c, posterior_cov_c) 
            
          }
          
          # split locations 
          mus_a_samp <- mus_samp[,1:M_a_samp,  drop = FALSE]
          if(M_na_samp != 0){
            mus_na_samp <- mus_samp[,(M_a_samp + 1):M_samp,  drop = FALSE]
          }
          
        }
        
      }
    
    }  
    
    if(update_weight){
      
      # 3. Update the allocated weights from a standard or Selberg Dirichlet prior (double check)
      filled_comps <- as.integer(names(table(allocs_samp)))
      alpha_post <- rep(alpha_prior, M_a_samp)
      alpha_post[filled_comps] <-  alpha_post[filled_comps] + table(allocs_samp )
      alpha_post <- c(alpha_post, rep(alpha_prior, M_na_samp))
      
      if(M_samp != 1){
        if(weights_prior == "standard"){
          weights_samp <- c(brms::rdirichlet(1, alpha_post))
        } else if(weights_prior == "selberg"){
          
         
          
          # fix for when values are so small distance is 0
          if(any(rdist(c(weights_samp * 1e100)) == 0)){
            log_demon <-  0
          }else{
            log_demon <-  sum(2 * gamma_samp * log(c(rdist(c(weights_samp * 1e100))) / 1e100))
          }
          
          w_prop <- c(brms::rdirichlet(1, alpha_post))

          if(any(rdist(c(w_prop * 1e100)) == 0)){
            log_numer <- 0
          }else{
            log_numer <- sum(2 * gamma_samp * log(c(rdist(c(w_prop * 1e100))) / 1e100))
          }
          
          # scaling_factor = 3
          # w_prop <- c(brms::rdirichlet(1, scaling_factor*alpha_post))
          # if(any(rdist(c(weights_samp * 1e100)) == 0)){
          #   log_demon <-  0
          # }else{
          #   log_demon <-  sum(2 * gamma_samp * log(c(rdist(c(weights_samp * 1e100))) / 1e100)) + sum(alpha_post * weights_samp) + sum(scaling_factor*alpha_post * w_prop)
          # }
          # 
          # 
          # if(any(rdist(c(w_prop * 1e100)) == 0)){
          #   log_numer <- 0
          # }else{
          #   log_numer <- sum(2 * gamma_samp * log(c(rdist(c(w_prop * 1e100))) / 1e100)) + sum(alpha_post * w_prop) + sum(scaling_factor*alpha_post * weights_samp)
          # }
          
          log_acrate_weights <- min(log_numer - log_demon, 1)
          if(log(runif(1)) < log_acrate_weights){
            weights_samp <-  w_prop
          }
        }
      }
    }
    
    # Update repulsion parameters
    if((g_to_z & update_gamma) | (g_time_z & update_gamma)){
      
      gamma_prop <- exp(rnorm(1, mean = log(gamma_samp) , sd = tuning_sd_gamma))
      
      log_acrate <-
        log_D(M_samp - 1, alpha_prior, beta_prior, gamma_samp) -
        log_D(M_samp - 1, alpha_prior, beta_prior, gamma_prop)
      if(g_to_z){
        log_acrate <-  log_acrate +
          # Only for 2d now:  
          # normalizing constants
          2 * log_G(M_samp, zetas_samp[1]) -
          2 * log_G(M_samp, gamma_prop * g_to_z_ratio) +
          # pairwise differences and prior and part coming from zeta likelihood, now for both zetas
          (-sum(g_to_z_ratio * mus_samp[1,]^2 / 2) + g_to_z_ratio * sum(log(c(rdist(c(mus_samp[1,]))))) -
             sum(g_to_z_ratio * mus_samp[2,]^2 / 2) + g_to_z_ratio * sum(log(c(rdist(c(mus_samp[2,])))))
           + 2 * sum(log(c(rdist(c(weights_samp))))) - rep_prior_rate) * (gamma_prop - gamma_samp) +
          rep_prior_shape * (log(gamma_prop) - log(gamma_samp))
      }
      
      log_acrate_gamma <- min(1, log_acrate)
      
      if(log(runif(1)) < log_acrate_gamma){
        gamma_samp <-  gamma_prop
      }
      
      if(g_to_z){
        zetas_samp <- rep(gamma_samp * g_to_z_ratio, D)
      }
      
      # else if(g_time_z){
      #   zetas_samp <- g_time_z_ratio / gamma_prop
      # }
      
    }else if((g_to_z|g_time_z) & !update_gamma){
      
      # Zeta cant be 0
      if(gamma_samp == 0){
        zetas_samp <- rep(1e-5, D)
      }else{
        if(g_to_z){
          zetas_samp <- rep(gamma_samp * g_to_z_ratio, D)
        }
        # else if(g_time_z){
        #   zetas_samp <- g_time_z_ratio / gamma_prop
        # }
      }
      
    }else{
      if(update_gamma){
        gamma_prop <- exp(rnorm(1, mean = log(gamma_samp) , sd = tuning_sd_gamma))
        log_acrate <- 
          # normalizing constants
          log_D(M_samp - 1, alpha_prior, beta_prior, gamma_samp) -
          log_D(M_samp - 1, alpha_prior, beta_prior, gamma_prop) +
          # pairwise differences and prior
          (2 * sum(log(c(rdist(c(weights_samp))))) - rep_prior_rate) * (gamma_prop - gamma_samp) +
          rep_prior_shape * (log(gamma_prop) - log(gamma_samp))
        log_acrate_gamma <- min(1, log_acrate)
        if(log(runif(1)) < log_acrate_gamma){
          gamma_samp <-  gamma_prop
        }
      }
      
      if(update_zetas){
        
        for(d in 1:D){
          zeta_prop <- exp(rnorm(1, mean = log(zetas_samp[d]) , sd = tuning_sd_zeta))
          log_acrate <-  
            # normalizing constants
            log_G(M_samp, zetas_samp[d]) -
            log_G(M_samp, zeta_prop) +
            # pairwise differences
            (-sum(mus_samp[d, ]^2 / 2) + sum(log(c(rdist(c(mus_samp[d, ]))))) - rep_prior_rate) * (zeta_prop - zetas_samp[d]) +
            rep_prior_shape * (log(zeta_prop) - log(zetas_samp[d]))
          log_acrate <- min(1, log_acrate)
          if(log(runif(1)) < log_acrate){
            zetas_samp[d] <-  zeta_prop
          }
        }
        
      } # update zeta
    } # g_to_z
    
    
    if(update_allocs & update_M ){
      # Birth-Death move
      # We define 1 to be a Birth-move and 0 a Death-move 
      # If there are no non-allocated components, we automatically have a Birth move
      birth_death <-ifelse(M_na_samp == 0, 1, sample(c(0,1), 1, prob = c(1 - p_b, p_b)) )
      #print(paste0(birth_death))
      actual_birth <- 0
      actual_death <- 0
      log_accept_birth <- 0
      log_accept_death <- 0
      # reset acceptance probabilities#
      # birth
      log_accept_birth_M <- 0
      log_accept_birth_weights <- 0
      log_accept_birth_locs <- 0
      log_accept_birth_bd_props <- 0
      log_accept_birth_bd_probs <- 0
      #death
      log_accept_death_M <- 0
      log_accept_death_weights <- 0
      log_accept_death_locs <- 0
      log_accept_death_bd_props <- 0
      log_accept_death_bd_probs <- 0
      if(birth_death == 1){
        # Start by proposing a new set of weights
        alpha_new <- alpha_prior
        weights_samp_new <- c(brms::rdirichlet(1, c(alpha_post, alpha_new))) #alpha_prior
        # Then a new location
        if(D == 1){
  
          log_accept_birth_locs <- 0
          
          if(repulsive_locations){
            #mu_new <- rnorm(1, mu_prop_means, sd = sqrt(mu_prop_Sigma))
            mu_new <- rnorm(1, 0, sd = sqrt(1 / zetas_samp))
            log_new_pwd_mu <-  sum(log(abs(cdist(mu_new, mus_samp))))
            log_accept_birth_locs <-
              log_accept_birth_locs +
              log_G(M_samp, zetas_samp) -
              log_G(M_samp + 1, zetas_samp) +
              #kernel ca
              #(zetas_samp / 2) * mu_new^2 +
              # pwd
              zetas_samp * (log_new_pwd_mu)
          }
          
        }else{
          #mu_new <- rmvnorm(1, mu_prop_means, sigma = mu_prop_Sigma)
          for(d in 1:D){
            mu_new <- rmvnorm(1, rep(0,D), sigma = diag(rep(1 / zetas_samp[d], D)))
            log_new_pwd_mu <-  sum(log(abs(cdist(mu_new[d], mus_samp[d,]))))
            
            log_accept_birth_locs <- 
              log_accept_birth_locs + 
              log_G(M_samp, zetas_samp[d]) -
              log_G(M_samp + 1, zetas_samp[d]) +
              #kernel cancels if variance is 1 / zetas, watch out for +/-
              #(zetas_samp[d] / 2) * mu_new[d]^2 +
              # pwd
              zetas_samp[d] * (log_new_pwd_mu)    
            
          }
        }
        # Also propose new Sigma
        if(D == 1){
          Sigma_new <- 1 / rgamma(1, shape = sigma2_prior_shape, rate = sigma2_prior_rate)
        }else{
          Sigma_new <- array(c(Sigmas_a_samp, Sigmas_na_samp), dim = c(D, D, M_samp))
        }
        # Pre-calculate pairwise differences for the new and old set of weights
        log_pwd_weights <- sum(log(c(rdist(c(weights_samp)))))
        log_new_pwd_weights <- sum(log(c(rdist(c(weights_samp_new)))))
        # Prior on number of components
        log_accept_birth_M <- log(M_prior) - log(M_samp) 
        # SDir prior on weights
        log_accept_birth_weights <- 
          # Note that the normalization constants take M - 1 as input!--> consider changing this
          log_D(M_samp - 1, alpha_prior, beta_prior, gamma_samp) -
          log_D(M_samp + 1 - 1, alpha_prior, beta_prior, gamma_samp) +
          # pwd
          2 * gamma_samp * (log_new_pwd_weights - log_pwd_weights)
        log_accept_birth_bd_props <- 
          lgamma(alpha_new) +
          lgamma(sum(c(alpha_post))) -  
          lgamma(sum(c(alpha_post, alpha_new))) -
          log(M_na_samp + 1) +
          #new part from the location proposal
          if(D == 1){
            if(repulsive_locations){
              # different if var is 1 / zetas
              # -0.5 * log(2 * pi * 1 / zetas_samp) + (1/ (2 * mu_prop_Sigma)) * (mu_new - mu_prop_means)^2
               0.5 * (log(2 * pi) - log(zetas_samp)) 
            }else{
              0
            }
          }else{
            if(repulsive_locations){
              #sum(log(dmvnorm_arma_mv(mu_new, mean = mu_prop_means, sigma = mu_prop_Sigma))) 
              # different if var is 1 / zetas
              sum(0.5 * (log(2 * pi) - log(zetas_samp))) 
            }else{
              0
            }
          }
          log_accept_birth_bd_probs <- ifelse(M_na_samp == 0, log(1 - p_b) - log(1), log(1 - p_b) - log(p_b) ) #log(1) = 0 so in case o Mna = 0 --> log(1 - p_b)
          #log_accept_birth_bd_probs <- log(1 - p_b) - log(p_b)
        # Adding all up
        log_accept_birth <- log_accept_birth_M + log_accept_birth_weights + log_accept_birth_locs + log_accept_birth_bd_props +  log_accept_birth_bd_probs
        # accept-reject 
        if(log(runif(1)) < log_accept_birth){
          actual_birth <- 1
          if(M_na_samp > 0){
            if(D == 1){
              mus_na_samp <- c(mus_na_samp, mu_new)
            }else{
              mus_na_samp <- cbind(mus_na_samp, t(mu_new)) 
            }
          }else{
            if(D == 1){
              mus_na_samp <- mu_new
            }else{
              mus_na_samp <- t(mu_new)
            }
          }
          M_na_samp <- M_na_samp + 1
          weights_samp <- weights_samp_new
          if(D == 1){
            Sigmas_na_samp <- c(Sigmas_na_samp, Sigma_new)
          }else{
            Sigmas_na_samp <- array(c(Sigmas_na_samp, Sigma_new), dim=c(D, D, M_na_samp) )
          }
          n_curr <- c(n_curr, 0)
          bd <- "b_b"
        }else{
          # indicator no birth happened
          actual_birth <- -1
          bd <- "b_no_b"
        }
      }else{ #Death-move
        # Start by choosing the component which should be closed
        na_close <- sample(M_na_samp, 1)
        weights_samp_new <- c(brms::rdirichlet(1, alpha_post[1:(M_samp - 1)]))
        # Pre-calculate pairwise differences of the component chosen to be closed and all others
        log_pwd <- sum(log(c(rdist(c(weights_samp)))))
        log_new_pwd <- sum(log(c(rdist(c(weights_samp_new)))))
        
        if(D == 1){
          
          log_accept_death_locs <- 0
          if(!repulsive_locations){
            log_new_pwd_mu <-  sum(log(abs(cdist(mus_samp[M_a_samp + na_close], mus_samp[-(M_a_samp + na_close)] ))))
            log_accept_death_locs <-
              log_accept_death_locs +
              # normalization constants
              log_G(M_samp, zetas_samp) -
              log_G(M_samp - 1, zetas_samp) -
              #kernel cancels if variance is 1 / zetas, watch out for +/-
              #(zetas_samp / 2) * mus_samp[d, M_a_samp + na_close]^2 -
              # pwd
              zetas_samp * (log_new_pwd_mu)
          }
        }else{
          
          for(d in 1:D){
            
            log_new_pwd_mu <-  sum(log(abs(cdist(mus_samp[d,M_a_samp + na_close], mus_samp[d, -M_a_samp + na_close]))))
            log_accept_death_locs <- 
              log_accept_death_locs +
              # normalization constants
              log_G(M_samp, zetas_samp[d]) -
              log_G(M_samp - 1, zetas_samp[d]) -
              #kernel cancels if variance is 1 / zetas, watch out for +/-
              #(zetas_samp[d] / 2) * mus_samp[d, M_a_samp + na_close]^2 -
              # pwd
              zetas_samp[d] * (log_new_pwd_mu)
            
          }
          
        }
        
        # Prior on number of components
        log_accept_death_M <- log(M_prior) + log(M_samp - 1) 
        # SDir prior on weights
        log_accept_death_weights <- 
          log_D(M_samp - 1, alpha_prior, beta_prior, gamma_samp) -
          log_D(M_samp - 1 - 1, alpha_prior, beta_prior, gamma_samp) +
          # pwd
          2 * gamma_samp * (log_new_pwd - log_pwd)
        # ratio of death & birth proposal densities
        log_accept_death_bd_props <-  
          log(M_na_samp) +
          lgamma(sum(c(alpha_post))) -
          lgamma(alpha_post[M_samp]) - 
          lgamma(sum(alpha_post[1:(M_samp - 1)])) +
          if(D == 1){
            if(repulsive_locations){
              # new part from the location proposal
              #-0.5 * log(2 * pi * 1 / zeta_samp) - (1/ (2 * mu_prop_var)) * (mus_samp[M_a_samp + na_close] - mu_prop_means)
              0.5 * (log(zetas_samp) - log(2 * pi))
            }else{
              0
            }
          }else{
            if(repulsive_locations){
              # new part from the location proposal
              #sum(log(dmvnorm_arma_mv(t(mus_samp[, M_a_samp + na_close, drop = FALSE]), mean = mu_prop_means, sigma = mu_prop_Sigma)))
              # different if var is 1 / zetas
              sum(0.5 * (log(zetas_samp) - log(1 / zetas_samp)))
            }else{
              0
            }
          }
        # ratio of birth and death probabilities
        log_accept_death_bd_probs <- - log(1 - p_b) + log(p_b)
        # Adding all up
        log_accept_death <- log_accept_death_M + log_accept_death_weights + log_accept_death_locs + log_accept_death_bd_props + log_accept_death_bd_probs
        # accept-reject 
        if(log(runif(1)) < log_accept_death){
          actual_death <- 1
          M_na_samp <- M_na_samp - 1
          weights_samp <- weights_samp_new
          if(D == 1){
            mus_na_samp <- mus_na_samp[-na_close]
            Sigmas_na_samp <- Sigmas_na_samp[-na_close]
          }else{
            mus_na_samp <- mus_na_samp[,-na_close, drop = FALSE]
            Sigmas_na_samp <- Sigmas_na_samp[,,-na_close, drop = FALSE] 
          }
          n_curr <- n_curr[-(M_a_samp + na_close)]
          bd <- "d_d"
        }else{
          actual_death <- -1
          bd <- "d_no_d"
        }
        
      }
      # Update total number of components
      M_samp <- M_a_samp + M_na_samp
      if(D == 1){
        mus_samp <- c(mus_a_samp, mus_na_samp)
        Sigmas_samp <- c(Sigmas_a_samp, Sigmas_na_samp)
      }else{
        mus_samp <- cbind(mus_a_samp, mus_na_samp)
        Sigmas_samp <- array(c(Sigmas_a_samp, Sigmas_na_samp), dim = c(D, D, M_samp))
      }
      
    }
    
    if(D == 1){
      if(length(mus_samp) != M_samp)stop("error after allocs")
    }else{
      if(ncol(mus_samp) != M_samp)stop("error after allocs")
    }
    
    if(update_post_dens){
      if(D == 1){
        post_dens <- rep(NA, length(y_line))
        for(i in 1:length(y_line)){
          
          post_dens[i] <- sum(weights_samp * dnorm(y_line[i],  mus_samp , sd = sqrt(Sigmas_samp )))
          
        }
      }else if(D == 2){
        
        post_dens <- post_density(y_grid,  weights_samp, mus_samp, Sigmas_samp)
      
        }
      
    }else{
      if(D == 1){
        post_dens <- rep(NA, length(y_line))
      }else if(D == 2){
        post_dens <- rep(NA, nrow(y_grid))
      }else{
        post_dens <- NA
      }
    }
    
    # Compute minimum pairwise differences and rankings
    
    #w_rank <- matrix(NA, 1, 4)
    #colnames(w_rank) <- c("correct", "avg", "min_pwd", "kendall")
    
    # w_rank[1,1] <- all(order(weights_samp[order(mus_samp)], decreasing = TRUE) == order(true_allocs, decreasing = TRUE))
    # w_rank[1,2] <- mean(order(weights_samp[order(mus_samp)], decreasing = TRUE) == order(true_allocs, decreasing = TRUE))
    # w_rank[1,3] <- min(c(rdist(weights_samp)))
    # w_rank[1,4] <- DistancePair(order(weights_samp), order(true_allocs, decreasing = TRUE))
     
    # w_rank[1,1] <- all(order(weights_samp) == order(true_weights))
    # w_rank[1,2] <- mean(order(weights_samp) == order(true_weights))
    # w_rank[1,3] <- min(c(rdist(weights_samp)))
    # w_rank[1,4] <- DistancePair(order(weights_samp), order(true_weights))
    
    #Save draws
    if (iter > n_burn && ((iter - n_burn) %% n_thin) == 0) {
      iter_aux <- (iter - n_burn) / n_thin
      M_a_out[iter_aux] <- M_a_samp
      M_na_out[iter_aux] <- M_na_samp
      alloc_out[iter_aux, ] <- allocs_samp
      #w_rank_out[iter_aux, ] <- w_rank
      
      if(D == 1){
        mu_out[iter_aux,] <- mus_samp#c(mus_a_samp, mus_na_samp)
        Sigma_out[iter_aux,] <- Sigmas_samp
        weights_out[iter_aux,] <- weights_samp
      }else{
        mu_out[[iter_aux]] <- cbind(mus_a_samp, mus_na_samp)
        Sigma_out[[iter_aux]] <- Sigmas_samp
        weights_out[[iter_aux]] <- weights_samp
      }
      
      gamma_out[iter_aux] <- gamma_samp
      zetas_out[iter_aux, ] <- zetas_samp
      bd_out[iter_aux] <- bd
      log_prior_mu_out[iter_aux] <- log_prior_mu
      log_lik_mu_out[iter_aux] <- log_lik_mu
      log_accept_death_out[iter_aux] <- log_accept_death
      log_accept_birth_out[iter_aux] <- log_accept_birth
      # log_accept_birth
      log_accept_birth_M_out[iter_aux] <- log_accept_birth_M 
      log_accept_birth_weights_out[iter_aux] <- log_accept_birth_weights 
      log_accept_birth_locs_out[iter_aux] <- log_accept_birth_locs 
      log_accept_birth_bd_props_out[iter_aux] <- log_accept_birth_bd_props
      
      post_dens_out[iter_aux,] <- post_dens
      # GE
      mu_norm_consts_out[iter_aux,] <- mu_norm_consts
      mu_kernels_out[iter_aux,] <- mu_kernels
      mu_pwd_terms_out[iter_aux,] <-  mu_pwd_terms
      #weights_order_out[iter_aux] <-  mu_pwd_terms
      
    }
    
    if(update_allocs){
      n_bd_out[iter,] <- c(birth_death, actual_birth, actual_death)
    }
    
    if (iter %% 100 == 0) {
      print(iter)
      #print(paste0(M_a_samp))
    }
  }
  
  results <- list(alloc_out = alloc_out, mu_out = mu_out, 
                  weights_out = weights_out, Sigma_out = Sigma_out,
                  M_a_out = M_a_out, M_na_out = M_na_out, 
                  n_bd_out = n_bd_out,
                  n_cs_out = n_cs_out,
                  gamma_out = gamma_out,  
                  zetas_out = zetas_out,
                  bd_out = bd_out, log_accept_death_out = log_accept_death_out,
                  log_accept_birth_out = log_accept_birth_out,
                  log_accept_birth_M_out = log_accept_birth_M_out,
                  log_accept_birth_weights_out = log_accept_birth_weights_out,
                  log_accept_birth_locs_out = log_accept_birth_locs_out,
                  log_accept_birth_bd_props_out= log_accept_birth_bd_props_out,
                  post_dens_out = post_dens_out,
                  mu_norm_consts_out = mu_norm_consts_out,
                  mu_kernels_out = mu_kernels_out,
                  mu_pwd_terms_out = mu_pwd_terms_out
                 # w_rank_out = w_rank_out
  )
  
  return(results)
}
