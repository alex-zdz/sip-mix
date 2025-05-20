# Define a new function for the birth-death move
birth_death_move <- function(update_allocs, update_M, M_na_samp, p_b, alpha_prior, alpha_post, D, 
                             repulsive_locations, mu_prop_means, mu_prop_Sigma, zetas_samp, M_samp, 
                             mus_samp, sigma2_prior_shape, sigma2_prior_rate, Sigmas_a_samp, 
                             Sigmas_na_samp, M_prior, gamma_samp, n_curr, M_a_samp) {
  if(update_allocs & update_M) {
    # Birth-Death move
    # We define 1 to be a Birth-move and 0 a Death-move 
    # If there are no non-allocated components, we automatically have a Birth move
    birth_death <- ifelse(M_na_samp == 0, 1, sample(c(0,1), 1, prob = c(1 - p_b, p_b)) )
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
          mu_new <- rnorm(1, 0, sd = sqrt(1 / zetas_samp))
          log_new_pwd_mu <-  sum(log(abs(cdist(mu_new, mus_samp))))
          log_accept_birth_locs <-
            log_accept_birth_locs +
            log_G(M_samp, zetas_samp) -
            log_G(M_samp + 1, zetas_samp) +
            zetas_samp * (log_new_pwd_mu)
        }
      } else {
        for(d in 1:D){
          mu_new <- rmvnorm(1, rep(0,D), sigma = diag(rep(1 / zetas_samp[d], D)))
          log_new_pwd_mu <-  sum(log(abs(cdist(mu_new[d], mus_samp[d,]))))
          log_accept_birth_locs <- 
            log_accept_birth_locs + 
            log_G(M_samp, zetas_samp[d]) -
            log_G(M_samp + 1, zetas_samp[d]) +
            zetas_samp[d] * (log_new_pwd_mu)    
        }
      }
      if(D == 1){
        Sigma_new <- 1 / rgamma(1, shape = sigma2_prior_shape, rate = sigma2_prior_rate)
      } else {
        Sigma_new <- array(c(Sigmas_a_samp, Sigmas_na_samp), dim = c(D, D, M_samp))
      }
      log_pwd_weights <- sum(log(c(rdist(c(weights_samp)))))
      log_new_pwd_weights <- sum(log(c(rdist(c(weights_samp_new)))))
      log_accept_birth_M <- log(M_prior) - log(M_samp) 
      log_accept_birth_weights <- 
        log_D(M_samp - 1, alpha_prior, beta_prior, gamma_samp) -
        log_D(M_samp + 1 - 1, alpha_prior, beta_prior, gamma_samp) +
        2 * gamma_samp * (log_new_pwd_weights - log_pwd_weights)
      log_accept_birth_bd_props <- 
        lgamma(alpha_new) +
        lgamma(sum(c(alpha_post))) -  
        lgamma(sum(c(alpha_post, alpha_new))) -
        log(M_na_samp + 1) +
        if(D == 1){
          if(repulsive_locations){
             0.5 * (log(2 * pi) - log(zetas_samp)) 
          } else {
            0
          }
        } else {
          if(repulsive_locations){
            sum(0.5 * (log(2 * pi) - log(zetas_samp))) 
          } else {
            0
          }
        }
      log_accept_birth_bd_probs <- ifelse(M_na_samp == 0, log(1 - p_b) - log(1), log(1 - p_b) - log(p_b) )
      log_accept_birth <- log_accept_birth_M + log_accept_birth_weights + log_accept_birth_locs + log_accept_birth_bd_props +  log_accept_birth_bd_probs
      if(log(runif(1)) < log_accept_birth){
        actual_birth <- 1
        if(M_na_samp > 0){
          if(D == 1){
            mus_na_samp <- c(mus_na_samp, mu_new)
          } else {
            mus_na_samp <- cbind(mus_na_samp, t(mu_new)) 
          }
        } else {
          if(D == 1){
            mus_na_samp <- mu_new
          } else {
            mus_na_samp <- t(mu_new)
          }
        }
        M_na_samp <- M_na_samp + 1
        weights_samp <- weights_samp_new
        if(D == 1){
          Sigmas_na_samp <- c(Sigmas_na_samp, Sigma_new)
        } else {
          Sigmas_na_samp <- array(c(Sigmas_na_samp, Sigma_new), dim=c(D, D, M_na_samp) )
        }
        n_curr <- c(n_curr, 0)
        bd <- "b_b"
      } else {
        actual_birth <- -1
        bd <- "b_no_b"
      }
    } else {
      na_close <- sample(M_na_samp, 1)
      weights_samp_new <- c(brms::rdirichlet(1, alpha_post[1:(M_samp - 1)]))
      log_pwd <- sum(log(c(rdist(c(weights_samp)))))
      log_new_pwd <- sum(log(c(rdist(c(weights_samp_new)))))
      if(D == 1){
        log_accept_death_locs <- 0
        if(!repulsive_locations){
          log_new_pwd_mu <-  sum(log(abs(cdist(mus_samp[M_a_samp + na_close], mus_samp[-(M_a_samp + na_close)] ))))
          log_accept_death_locs <-
            log_accept_death_locs +
            log_G(M_samp, zetas_samp) -
            log_G(M_samp - 1, zetas_samp) -
            zetas_samp * (log_new_pwd_mu)
        }
      } else {
        for(d in 1:D){
          log_new_pwd_mu <-  sum(log(abs(cdist(mus_samp[d,M_a_samp + na_close], mus_samp[d, -M_a_samp + na_close]))))
          log_accept_death_locs <- 
            log_accept_death_locs +
            log_G(M_samp, zetas_samp[d]) -
            log_G(M_samp - 1, zetas_samp[d]) -
            zetas_samp[d] * (log_new_pwd_mu)
        }
      }
      log_accept_death_M <- log(M_prior) + log(M_samp - 1) 
      log_accept_death_weights <- 
        log_D(M_samp - 1, alpha_prior, beta_prior, gamma_samp) -
        log_D(M_samp - 1 - 1, alpha_prior, beta_prior, gamma_samp) +
        2 * gamma_samp * (log_new_pwd - log_pwd)
      log_accept_death_bd_props <-  
        log(M_na_samp) +
        lgamma(sum(c(alpha_post))) -
        lgamma(alpha_post[M_samp]) - 
        lgamma(sum(alpha_post[1:(M_samp - 1)])) +
        if(D == 1){
          if(repulsive_locations){
            0.5 * (log(zetas_samp) - log(2 * pi))
          } else {
            0
          }
        } else {
          if(repulsive_locations){
            sum(0.5 * (log(zetas_samp) - log(2 * pi)))
          } else {
            0
          }
        }
    }
  }
  # Return the modified variables
  return(list(
    actual_birth = actual_birth,
    actual_death = actual_death,
    M_na_samp = M_na_samp,
    weights_samp = weights_samp,
    mus_na_samp = mus_na_samp,
    Sigmas_na_samp = Sigmas_na_samp,
    n_curr = n_curr,
    bd = bd
  ))
}