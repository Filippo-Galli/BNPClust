library(ggplot2)

# Likelihood components for 2 points definition ----

cohesion_term <- function(distance, delta_1, alpha, beta) {
  # Cohesion term calculation
  term1 <- distance^(delta_1 - 1) / gamma(delta_1)
  term2 <- (beta)^(alpha) / gamma(alpha)
  term3 <- gamma(delta_1 + alpha) / (distance + beta)^(delta_1 + alpha)

  term1 * term2 * term3
}

repulsion_term <- function(distance, delta_2, zeta, gamma) {
  # Repulsion term calculation
  term1 <- distance^(delta_2 - 1) / gamma(delta_2)
  term2 <- (gamma)^(zeta) / gamma(zeta)
  term3 <- gamma(delta_2 + zeta) / (distance + gamma)^(delta_2 + zeta)

  term1 * term2 * term3
}

distances <- seq(0.1, 2.5, by = 0.1)

# Cohesion Analysis ----
delta_1_values <- seq(0.1, 1, by = 0.2)
alpha_values <- seq(1, 10, by = 2)
beta_values <- seq(1, 10, by = 2)

## Plotting Cohesion Term for different delta_1 values for all distances

cohesion_data <- expand.grid(distance = distances, delta_1 = delta_1_values)
cohesion_data$cohesion_value <- mapply(cohesion_term, 
                    cohesion_data$distance, 
                    cohesion_data$delta_1, 
                    alpha = 1, beta = 1)

ggplot(cohesion_data, aes(x = distance, y = cohesion_value, color = factor(delta_1))) +
  geom_line() +
  labs(x = "Distance", 
     y = "Cohesion Term",
     title = "Cohesion Term vs Distance - alpha = 1, beta = 1",
     color = "delta_1") +
  theme_minimal()

### Save the plot
# ggsave("cohesion_delta1_plot.png", width = 8, height = 6)

## Plotting Cohesion Term for different alpha values for all distances

cohesion_data <- expand.grid(distance = distances, alpha = alpha_values)
cohesion_data$cohesion_value <- mapply(cohesion_term, 
                    cohesion_data$distance, 
                    0.9, 
                    alpha = cohesion_data$alpha, beta = 1)

ggplot(cohesion_data, aes(x = distance, y = cohesion_value, color = factor(alpha))) +
  geom_line() +
  labs(x = "Distance",
     y = "Cohesion Term",
     title = "Cohesion Term vs Distance - delta_1 = 0.9, beta = 1",
     color = "alpha") +
  theme_minimal()

# ggsave("cohesion_alpha_plot.png", width = 8, height = 6)

## Plotting Cohesion Term for different beta values for all distances
cohesion_data <- expand.grid(distance = distances, beta = beta_values)
cohesion_data$cohesion_value <- mapply(cohesion_term, 
                    cohesion_data$distance, 
                    0.9, 
                    alpha = 1, beta = cohesion_data$beta)

ggplot(cohesion_data, aes(x = distance, y = cohesion_value, color = factor(beta))) +
  geom_line() +
  labs(x = "Distance",
       y = "Cohesion Term",
       title = "Cohesion Term vs Distance - delta_1 = 0.9, alpha = 1",
       color = "beta") +
  theme_minimal()

ggsave("cohesion_beta_plot.png", width = 8, height = 6)

# Simultaneous increment of alpha/2 and beta
cohesion_data <- expand.grid(distance = distances, beta = beta_values)
cohesion_data$cohesion_value <- mapply(cohesion_term, 
                    cohesion_data$distance, 
                    0.9, 
                    alpha = cohesion_data$beta/2, beta = cohesion_data$beta)

ggplot(cohesion_data, aes(x = distance, y = cohesion_value, color = factor(beta))) +
  geom_line() +
  labs(x = "Distance",
       y = "Cohesion Term",
       title = "Cohesion Term vs Distance - delta_1 = 0.9 and alpha = beta/2",
       color = "beta") +
  theme_minimal()

ggsave("cohesion_alpha_beta_plot.png", width = 8, height = 6)

# Repulsion Analysis ----
distances <- seq(0.1, 20, by = 0.1)
delta_2_values <- seq(5, 21, by = 4)
zeta_values <- seq(1, 10, by = 2)
gamma_values <- seq(1, 10, by = 2)

## Plotting Repulsion Term for different delta_2 values for all distances
repulsion_data <- expand.grid(distance = distances, delta_2 = delta_2_values)
repulsion_data$repulsion_value <- mapply(repulsion_term,
                      repulsion_data$distance,
                      repulsion_data$delta_2,
                      zeta = 1, gamma = 1)

ggplot(repulsion_data, aes(x = distance, y = repulsion_value, color = factor(delta_2))) +
  geom_line() +
  labs(x = "Distance",
       y = "Repulsion Term",
       title = "Repulsion Term vs Distance - zeta = 1, gamma = 1",
       color = "delta_2") +
  theme_minimal()

ggsave("repulsion_delta2_plot.png", width = 8, height = 6)

## Plotting Repulsion Term for different zeta values for all distances
repulsion_data <- expand.grid(distance = distances, zeta = zeta_values)
repulsion_data$repulsion_value <- mapply(repulsion_term,
                      repulsion_data$distance,
                      10,
                      zeta = repulsion_data$zeta, gamma = 1)  
ggplot(repulsion_data, aes(x = distance, y = repulsion_value, color = factor(zeta))) +
  geom_line() +
  labs(x = "Distance",
        y = "Repulsion Term", 
        title = "Repulsion Term vs Distance - delta_2 = 10, gamma = 1",
        color = "zeta") +
  theme_minimal()

ggsave("repulsion_zeta_plot.png", width = 8, height = 6)

## Plotting Repulsion Term for different gamma values for all distances
repulsion_data <- expand.grid(distance = distances, gamma = gamma_values)
repulsion_data$repulsion_value <- mapply(repulsion_term,
                      repulsion_data$distance,
                      10,
                      zeta = 1, gamma = repulsion_data$gamma)
ggplot(repulsion_data, aes(x = distance, y = repulsion_value, color = factor(gamma))) +
  geom_line() +
  labs(x = "Distance",
       y = "Repulsion Term",
       title = "Repulsion Term vs Distance - delta_2 = 10, zeta = 1",
       color = "gamma") +
  theme_minimal()

ggsave("repulsion_gamma_plot.png", width = 8, height = 6)

## Simultaneous increment of zeta and gamma
repulsion_data <- expand.grid(distance = distances, gamma = gamma_values)
repulsion_data$repulsion_value <- mapply(repulsion_term,
                      repulsion_data$distance,
                      10,
                      zeta = repulsion_data$gamma*10, gamma = repulsion_data$gamma)
ggplot(repulsion_data, aes(x = distance, y = repulsion_value, color = factor(gamma))) +
  geom_line() +
  labs(x = "Distance",
       y = "Repulsion Term",
       title = "Repulsion Term vs Distance - delta_2 = 10 and zeta = gamma*10",
       color = "gamma") +
  theme_minimal()

ggsave("repulsion_zeta_gamma_plot_delta.png", width = 8, height = 6)

# Density of v ----
conditional_density_V <- function(v, K, n, a, sigma, tau) {
  
  # Pre-compute exp(v)
  exp_v <- exp(v)
  
  # Pre-compute frequently used values
  a_over_sigma <- a / sigma
  tau_power_sigma <- tau^sigma
  
  # Compute log density components
  # log(e^{vn}) = vn
  term1 <- v * n
  
  # log((e^v + τ)^{n-a|π|}) = -(n - a*K) * log(e^v + τ)
  term2 <- -(n - a * K) * log(exp_v + tau)
  
  # -(a/σ)((e^v+τ)^σ - τ^σ)
  term3 <- -a_over_sigma * ((exp_v + tau)^sigma - tau_power_sigma)
  
  # Return density
  exp(term1 + term2 + term3)
}

v_values <- seq(0, 20, length.out = 2000)
pi_values <- seq(10, 40, by = 5)
n = 100
a_values <- seq(0.1, 0.9, by = 0.2)
sigma_values <- seq(0.1, 0.9, by = 0.2)
tau_values <- seq(0.1, 2, by = 0.5)

## Plotting log density for different pi values
log_density_data <- expand.grid(v = v_values, a = a_values)
log_density_data$log_density_value <- mapply(conditional_density_V,
                                              log_density_data$v,
                                              7,
                                              100, 
                                              log_density_data$a, 0.7, 1)
ggplot(log_density_data, aes(x = v, y = log_density_value, color = factor(a))) +
  geom_line() +
  labs(x = "v",
       y = "Log Density of v",
       title = "Density of V unormalized - K = 7, n = 100, sigma = 0.7, tau = 1",
       color = "a values") +
  theme_minimal()

ggsave("log_density_a_plot.png", width = 8, height = 6)

## Plotting log density for different pi values
log_density_data <- expand.grid(v = v_values, pi = pi_values)
log_density_data$log_density_value <- mapply(conditional_density_V,
                                              log_density_data$v,
                                              log_density_data$pi,
                                              100, 0.01, 0.7, 1)
ggplot(log_density_data, aes(x = v, y = log_density_value, color = factor(pi))) +
  geom_line() +
  labs(x = "v",
       y = "Log Density of v",
       title = "Density of V unormalized - n = 100, a = 0.01, sigma = 0.7, tau = 1",
       color = "Number of clusters (|π|)") +
  theme_minimal()

ggsave("log_density_pi_plot.png", width = 8, height = 6)
