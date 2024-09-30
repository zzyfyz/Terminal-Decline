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
  cluster <- 50
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
  
  b <- 0.03
  c <- 0.05
  lambda0 <- 0.04
  gamma <- 1.2
  sigma_u <- 5
  sigma_b <- 8
  sigma_e <- 3
  
  # Fixed effects covariates
  x1 <- rbinom(n, 1, 0.5)
  x2 <- scale(rnorm(n, mean = 70, sd = 10), center = TRUE, scale = FALSE)  # Center x2
  time_points <- seq(0, time, by = 3)
  subject_cluster <- rep(1:cluster, each = cluster_subj)
  treatment_clusters <- sample(1:cluster, size = cluster / 2, replace = FALSE)
  treatment <- ifelse(subject_cluster %in% treatment_clusters, 0, 1)
  
  # Random effects
  ui <- rnorm(cluster, mean = 0, sd = sigma_u)
  bi <- rnorm(n, mean = 0, sd = sigma_b)
  epsiloni <- mvrnorm(n, mu = rep(0, length(time_points)), Sigma = diag(sigma_e^2, length(time_points)))
  
  # Cox frailty model for survival data assuming Weibull distribution
  linear <- alpha11 * x1 + alpha12 * x2 + b * bi + c * ui[subject_cluster]
  lambda <- lambda0 * exp(linear)
  U <- runif(n)
  
  # Simulate survival times
  survival_times <- (-log(U) / lambda)^(1 / gamma)
  censoring_times <- time
  observed_times <- pmin(survival_times, censoring_times)
  status <- as.numeric(survival_times <= censoring_times)
  
  # Create a data frame for subjects and time points
  longitudinal_data <- expand.grid(subject = 1:n, time = time_points)
  
  # Merge subject-level data with longitudinal data
  longitudinal_data <- longitudinal_data %>%
    mutate(
      x1 = x1[subject],
      x2 = x2[subject],
      cluster = subject_cluster[subject],
      treatment = treatment[subject],
      survival_time = survival_times[subject],
      observed_time = observed_times[subject],
      status = status[subject],
      backward_time = survival_time - time,
      bi = bi[subject],
      ui = ui[cluster]
    )
  
  # Calculate longitudinal measurements using vectorized operations
  longitudinal_data <- longitudinal_data %>%
    mutate(
      measurement = ifelse(
        backward_time > 0,
        alpha00 + x1 * alpha01 + x2 * alpha02 + backward_time * alpha03 +
          treatment * alpha04 + bi + ui + epsiloni[cbind(subject, match(time, time_points))],
        NA
      )
    )
  
  # Convert to wide format
  longitudinal_data_wide <- longitudinal_data %>%
    select(subject, time, measurement) %>%
    pivot_wider(names_from = time, values_from = measurement, names_prefix = "time_")
  
  # Merge subject-level data with the wide-format longitudinal data
  final_data <- longitudinal_data %>%
    distinct(subject, .keep_all = TRUE) %>%
    select(subject, x1, x2, cluster, treatment, survival_time, observed_time, status) %>%
    left_join(longitudinal_data_wide, by = "subject")
  
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

