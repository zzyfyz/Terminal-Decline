library(MASS)
library(survival)
library(tidyr)
library(dplyr)
dirg <- "C:/Users/feiyi/OneDrive/Desktop/Katie/Terminal-Decline"
setwd(dirg)
num_simulations <- 100
# Simulation loop
for (sim in 1:num_simulations) {
  set.seed(123 + sim)  
  # Parameters
  cluster <- 40
  cluster_subj <- 20
  n <- cluster * cluster_subj
  time <- 12
  alpha00 <- 10
  alpha01 <- 5
  alpha02 <- 0.3
  alpha03 <- 2 
  alpha04 <- 10 
  alpha11 <- -1
  alpha12 <- 0.05 
  #alpha13 <- 2
  
  b <- 1
  c <- 1
  lambda0 <- 0.08
  gamma <- 1.1
  sigma_u <- 1
  sigma_b <- 1
  sigma_e <- 1
  
  # Fixed effects covariates
  x1 <- rbinom(n, 1, 0.5)
  x2 <- rnorm(n, mean = 70, sd = 10)
  x2 <- x2-mean(x2)
  time_points <- seq(0,time, by=3)
  subject_cluster <- rep(1:cluster, each = cluster_subj)
  treatment_clusters <- sample(1:50, size = 25, replace = FALSE)
  treatment <- ifelse(subject_cluster %in% treatment_clusters, 0, 1)
  
  # Random effects
  ui <- rnorm(cluster, mean = 0, sd = sigma_u)
  bi <- rnorm(n, mean = 0, sd = sigma_b)
  epsiloni <- mvrnorm(n, mu = rep(0, time), Sigma = diag(sigma_e*sigma_e, time))
  
  # Cox frailty model for survival data assuming Weibull distribution
  linear <- alpha11 * x1 + alpha12 * x2  + b * bi + c * ui[subject_cluster]
  lambda <- lambda0 * exp(linear)
  U <- runif(n)
  
  # Simulate survival times
  survival_times <- (-log(U)/lambda)^(1/gamma)
  censoring_times <- time 
  observed_times <- pmin(survival_times, censoring_times)
  status <- as.numeric(survival_times <= censoring_times)
  
  # Mixed model for longitudinal data, modeling backward from death
  longitudinal_data <- data.frame()
  
  for (i in 1:n) {
    for (ind in seq_along(time_points)) {
      t = time_points[ind]
      backward_time <- survival_times[i] - t
      if (backward_time > 0) {
        measurement <- alpha00 + x1[i] * alpha01 + x2[i] * alpha02 + backward_time * alpha03 + treatment[i] * alpha04 + bi[i] + ui[subject_cluster[i]] + epsiloni[i, ind]
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
  mask <- !is.na(as.matrix(longitudinal_data_wide[, grep("time_", colnames(longitudinal_data_wide))]))
  
  longitudinal_data_wide[is.na(longitudinal_data_wide)] <- 0
  
  final_data <- longitudinal_data_wide %>%
    left_join(data.frame(subject = 1:n, x1, x2, cluster = subject_cluster, treatment = treatment, time = observed_times, death = survival_times, status),
              by = "subject")
  
  mask_df <- as.data.frame(mask)
  filename <- paste0("mask.", sim-1, ".csv")
  write.csv(mask_df, file = filename, row.names = FALSE)
  
  filename <- paste0("sim.data.", sim-1, ".csv")
  write.csv(final_data, file = filename, row.names = FALSE)
  
}

