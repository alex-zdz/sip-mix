library(ggplot2)
library(dplyr)

# Define the derivative function
derivative_variance <- function(gamma, M, alpha) {
  num1 <- (M - 2) * (M - 1) * alpha
  inner1 <- (M - 2) * (M - 1) * gamma
  inner2 <- 2 * (M - 2) * (M - 1) * gamma + (4 * M - 3) * alpha + 1
  inner3 <- alpha * (M * ((2 * M - 3) * alpha + 1) - 2)
  numerator <- -num1 * (inner1 * inner2 + inner3)
  
  denom1 <- ((M - 2) * (M - 1) * gamma + M * alpha)^3
  denom2 <- ((M - 2) * (M - 1) * gamma + M * alpha + 1)^2
  denominator <- denom1 * denom2
  
  return(numerator / denominator)
}

# Set values
gamma_vals <- seq(0.01, 5, length.out = 500)
M_vals <- c(3, 5, 7, 10)  # Add/remove M values as needed
alphas <- c(0.1, 1, 3, 5, 10)

# Create a data frame with all combinations
df <- expand.grid(gamma = gamma_vals, M = M_vals, alpha = alphas) %>%
  mutate(deriv = mapply(derivative_variance, gamma, M, alpha = alpha))

# Plot using ggplot
ggplot(df, aes(x = gamma, y = deriv, col = factor(alpha))) +
  geom_line( size = 1) +
  facet_wrap(~ M, scales = "free_y") + #, ncol = 2
  labs(title = bquote("Derivative of Variance with respect to " ~ gamma ~ " for " ~ alpha == .(alpha)),
       x = expression(gamma),
       y = expression(frac(d, d*gamma)*"Var")) +
  theme_minimal(base_size = 14)


library(tidyr)

# 1. Repulsive Dirichlet variance
variance_repulsive <- function(alpha, M, gamma) {
  k <- alpha * M + (M - 1) * (M - 2) * gamma
  if (k <= 0 || k + 1 <= 0) return(NA)
  v <- (alpha / k) * ((1 - (alpha / k)) / (k + 1))
  return(v)
}

# 2. Standard Dirichlet variance via E[X^2] - E[X]^2
variance_standard_mean_method <- function(alpha, M) {
  if (M * alpha <= 0 || M * alpha + 1 <= 0) return(NA)
  mean_sq <- (alpha / (M * alpha))^2
  mean_of_sq <- ((alpha + 1) * alpha) / ((M * alpha)^2 * (M * alpha + 1))
  return(mean_of_sq - mean_sq)
}

# 3. Wikipedia formula for symmetric Dirichlet
variance_standard_wiki <- function(alpha, M) {
  alpha0 <- M * alpha
  if (alpha0 <= 0 || alpha0 + 1 <= 0) return(NA)
  return((alpha * (alpha0 - alpha)) / (alpha0^2 * (alpha0 + 1)))
}

# Grid of alpha and M values
alpha_vals <- seq(0.1, 20, length.out = 100)
M_vals <- c(3, 5, 7, 10)
gamma_vals <- c( 0.5, 1, 3, 5)

# Create full grid
grid <- expand.grid(alpha = alpha_vals, M = M_vals, gamma = gamma_vals)

df <- grid %>%
  mutate(
    repulsive = mapply(variance_repulsive, alpha, M, gamma),
    dirichlet = mapply(variance_standard_wiki, alpha, M),
  ) %>%
  pivot_longer(cols = c(repulsive, dirichlet), 
               names_to = "type",
               values_to = "variance")

# Plot
ggplot(df, aes(x = alpha, y = variance, color = type)) +
  geom_line(size = 1) +
  facet_grid(gamma ~ M, scales = "free_y",
             labeller = label_both) +
  labs(
    title = bquote("Variance of " ~ w[j] ~ " in Repulsive and Standard Dirichlet Models"),
    x = expression(alpha),
    y = expression(Var(w[j])),
    color = "Type"
  ) +
  scale_color_npg() +
  theme_minimal(base_size = 14)


df %>%filter(gamma != 1, M != 7) %>%
ggplot( aes(x = alpha, y = variance, color = type)) +
  geom_line(size = 1) +
  facet_grid(M ~ gamma, scales = "free_y",
             labeller = label_both) +
  labs(
    title = bquote("Variance of " ~ w[j] ~ " in Repulsive and Standard Dirichlet Models"),
    x = expression(alpha),
    y = expression(Var(w[j])),
    color = "Type"
  ) +
  scale_color_npg() +
  theme_minimal(base_size = 14)
