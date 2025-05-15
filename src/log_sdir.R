log_sdir <- function(x, alpha, beta, gamma){
  p <- length(x)
  if(p != 1){
    log_norm_const <- log_D(p, alpha, beta, gamma)
    # Pairwise differences
    log_pwd <- sum(log(c(rdist(x))))
    log_kernel <- (alpha - 1) * sum(log(x)) + (beta - 1) * log(1 - sum(x)) + 2 * gamma * log_pwd
    return( - log_norm_const + log_kernel) 
  }else{
    return(0)
  }
}

log_D <- function(p, alpha, beta, gamma){
  ind <- 1:p
  if(p != 0){
    log_D <- lgamma(beta) - lgamma(beta + alpha * p + gamma * p * (p - 1)) + 
      sum( lgamma(1 + ind * gamma) + lgamma(alpha + gamma *(p - ind)) - lgamma(1 + gamma) )
  }else{
    log_D <- 0
  }    
  return(log_D)
}  

# Gaussian ensemble

log_G <- function(p, zeta){
  #if(zeta != 0){
    ind <- 0:(p-1)
    log_G <- (- p / 2 - zeta * p *(p - 1) / 4 ) * log(zeta) + p / 2 * log(2 * pi) + 
      sum(lgamma(1 + (ind + 1) * zeta / 2 ) - lgamma(1 + zeta / 2 ) ) 
  #}else{
  #  log_G <- 0
  #}
  return(log_G)
}

G <- function(p, zeta){
  
  ind <- 0:(p-1)
  G <- zeta^(-p / 2 - zeta * p * (p - 1) / 4) * (2 * pi) ^ (p / 2) *
  prod( gamma(1 + (1 + ind) * (zeta / 2)) / gamma(1 + zeta / 2) )  
  
  return(G)
  
}

ge_dens <- function(x, z){
  
  #log_dens <- - (z / 2) * x[1]^2 + z * sum(log(cdist(x[1],x[2])))
  p <- length(x)
  dens <- (1 / G(p = 2, zeta = z)) * prod(exp(- (z / 2) * x^2 ))  * prod(c(rdist(x)))^z
  return(dens)
  
}



sdir <- function(y, alpha, beta, gamma){
  
  p <- length(y)
  abs(prod(c(rdist(y))))^(2*gamma) * (1 - sum(y))^(beta - 1) * prod(y^(alpha - 1)) *
    1 / sdir_norm(alpha, beta, gamma, p)
  
}



sdir_2 <- function(y, alpha, gamma){
  
  p <- length(y)
  abs(prod(c(rdist(y))))^(2*gamma) * prod(y^(alpha - 1)) *
    1 / sdir_norm(alpha, beta = alpha, gamma, p - 1)
  
}

log_sdir_2 <- function(y, alpha, gamma){
  
  p <- length(y)
  log_norm_const <- log_D(p - 1, alpha, beta = alpha, gamma)
  log_pwd <- sum(log(c(rdist(y))))
  log_kernel <- (alpha - 1) * sum(log(y)) + 2 * gamma * log_pwd
  return( - log_norm_const + log_kernel)
  
}


sdir_norm <- function(alpha, beta, gamma, p){
  
  ind_p <- c(1:p)
  first <- gamma(beta) / gamma(beta + alpha * p + gamma * p * (p - 1))
  second <- prod(gamma(alpha + gamma * (p - ind_p)) * gamma(1 + gamma * ind_p) / gamma(1 + gamma))
  return(first * second)
  
}

dir_norm <- function(alphas){
  
  dir_norm <- gamma(sum(alphas)) / prod(gamma(alphas))
  return(dir_norm)
  
}


log_dir_norm <- function(alphas){
  
  dir_norm <- lgamma(sum(alphas)) -  sum(lgamma(alphas))
  
  return(dir_norm)
  
}

# Selberg Beta distribution

sbe <- function(x, alpha, beta, zeta){
  
  M <- length(x)
  kernel <- prod(x^(alpha - 1) * (1 - x)^(beta - 1))
  pwd <- abs(prod(c(rdist(x))))^(2*zeta)
  B <- sbe_norm(M, alpha, beta, zeta)
  return((1 / B) * kernel * pwd)
  
}

sbe_norm <- function(M, alpha, beta, zeta){
  
  ind <- 1:M
  numer <- gamma(alpha + (ind - 1) * zeta) * gamma(beta + (ind - 1) * zeta) * gamma(1 + ind * zeta)  
  demon <- gamma(alpha + beta + (M + ind - 2) * zeta) * gamma(1 + zeta)
  return(prod(numer / demon))
  
}

# Log-domain

log_sbe <- function(x, alpha, beta, zeta) {
  M <- length(x)
  log_kernel <- sum((alpha - 1) * log(x) + (beta - 1) * log(1 - x))
  log_pwd <- 2 * zeta * sum(log(c(rdist(x))))
  log_B <- log_sbe_norm(M, alpha, beta, zeta)
  return(log_kernel + log_pwd - log_B)
}


log_sbe_norm <- function(M, alpha, beta, zeta) {
  ind <- 1:M
  log_numer <- lgamma(alpha + (ind - 1) * zeta) + lgamma(beta + (ind - 1) * zeta) + lgamma(1 + ind * zeta)
  log_demon <- lgamma(alpha + beta + (M + ind - 2) * zeta) + lgamma(1 + zeta)
  log_density <- sum(log_numer - log_demon)
  return(log_density)
}







