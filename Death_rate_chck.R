library(MASS)
library(survival)
library(tidyr)
library(dplyr)
dirg <- "C:/Users/feiyi/OneDrive/Desktop/Katie/Terminal-Decline/Spline"
setwd(dirg)
num_simulations <- 100
# Simulation loop
for (sim in 1:num_simulations) {
  set.seed(123 + sim)  
  # Parameters
  cluster <- 20
  cluster_subj <- 50
  n <- cluster * cluster_subj
  time <- 24
  alpha00 <- 30
  alpha01 <- 1
  alpha02 <- 0.9
  alpha03 <- -30
  alpha04 <- 0.2 
  alpha05 <- 30 
  alpha06 <- -0.23 
  alpha07 <- -0.92 
  
  alpha11 <- 0.2
  alpha12 <- -0.01
  #alpha13 <- 2
  
  b <- 0.03
  c <- 0.02
  lambda0 <- 0.05
  gamma <- 2.2
  sigma_u <- 5
  sigma_b <- 6
  sigma_e <- 4
  
  # Fixed effects covariates
  x1 <- rbinom(n, 1, 0.5)
  x2 <- runif(n, 100, 150)
  time_points <- seq(1,time, by=6)
  subject_cluster <- rep(1:cluster, each = cluster_subj)
  treatment_clusters <- sample(1:cluster, size = cluster/2, replace = FALSE)
  treatment <- ifelse(subject_cluster %in% treatment_clusters, 0, 1)
  
  # Random effects
  ui <- rnorm(cluster, mean = 0, sd = sigma_u)
  bi <- rnorm(n, mean = 0, sd = sigma_b)
  epsiloni <- mvrnorm(n, mu = rep(0, time), Sigma = diag(sigma_e*sigma_e, time))
  
  # Cox frailty model for survival data assuming Weibull distribution
  linear <- alpha11 * x1 + alpha12 * x2 + b * bi + c * ui[subject_cluster]
  lambda <- lambda0 * exp(linear)
  U <- runif(n)
  
  # Simulate survival times
  survival_times <- (-log(U) / lambda)^(1 / gamma)
  censoring_times <- time
  observed_times <- pmin(survival_times, censoring_times)
  status <- as.numeric(survival_times <= censoring_times)
  
  # Mixed model for longitudinal data, modeling backward from death
  longitudinal_data <- data.frame()
  
  for (i in 1:n) {
    for (ind in seq_along(time_points)) {
      t <- time_points[ind]
      backward_time <- survival_times[i] - t
      if (backward_time > 0) {
        measurement <- alpha00 +
          (alpha03 / (1 + alpha04 * backward_time)) + 
          treatment[i] * alpha05 * exp(alpha06 * backward_time + alpha07) + 
          alpha01 * x1[i] + alpha02 * x2[i] + 
          bi[i] + ui[subject_cluster[i]] + epsiloni[i, ind]
      } else {
        measurement <- NA
      }
      longitudinal_data <- rbind(longitudinal_data, data.frame(subject = i, time = t, measurement = measurement))
    }
  }
  
  # Join the longitudinal and survival data
  longitudinal_data <- longitudinal_data %>%
    left_join(data.frame(subject = 1:n, observed_time = observed_times), by = "subject") %>%
    mutate(measurement = ifelse(time > observed_time, NA, measurement)) %>%
    select(-observed_time)
  
  # Create a mask for the missing data
  longitudinal_data_wide <- longitudinal_data %>%
    pivot_wider(names_from = time, values_from = measurement, names_prefix = "time_")
  
  # Add covariates and survival data to the final dataset without generating duplicates
  final_data <- longitudinal_data_wide %>%
    left_join(
      data.frame(subject = 1:n, x1, x2, cluster = subject_cluster, treatment = treatment, observed_time = observed_times, survival_time = survival_times, status = status),
      by = "subject"
    )
  
  # Create a mask for missing data
  measurement_columns <- grep("time_", names(longitudinal_data_wide), value = TRUE)
  mask <- !is.na(as.matrix(longitudinal_data_wide[, measurement_columns]))
  
  final_data[is.na(final_data)] <- 0
  
  mask_df <- as.data.frame(mask)
  filename <- paste0("mask.", sim-1, ".csv")
  write.csv(mask_df, file = filename, row.names = FALSE)
  
  filename <- paste0("sim.data.", sim-1, ".csv")
  write.csv(final_data, file = filename, row.names = FALSE)
  
}