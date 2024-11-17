library(dplyr)
library(rstan)
library(parallel)
library(splines2)
library(lme4)
library(ggplot2)

set.seed(999)  
# Parameters
cluster <- 50
cluster_subj <- 20
n <- cluster * cluster_subj
time <- 6
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
time_points <- seq(1,time, by=1)
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

# Create a mask to indicate missing values
mask <- !is.na(as.matrix(longitudinal_data_wide[, grep("time_", colnames(longitudinal_data_wide))]))

# Replace NA values with 0 in the longitudinal data
longitudinal_data_wide[is.na(longitudinal_data_wide)] <- 0

# Add covariates and survival data to the final dataset without generating duplicates
final_data <- longitudinal_data_wide %>%
  left_join(
    data.frame(subject = 1:n, x1, x2, cluster = subject_cluster, treatment = treatment, observed_time = observed_times, survival_time = survival_times, status = status),
    by = "subject"
  )

dat_death <- subset(final_data, status==1)
observed_times_death <- dat_death$observed_time

#Calculate backward times for each subject at each time point
backward_time_list <- lapply(observed_times_death, function(obs_time) {
  backward_times <- obs_time - time_points
  backward_times[backward_times < 0] <- NA  # Set negative backward times to NA
  return(backward_times)
})

#Flatten the list of backward times into a vector
backward_time_vector <- unlist(backward_time_list)

#Remove NA values from the vector
backward_time_vector <- backward_time_vector[!is.na(backward_time_vector)]

num_knots <- 8  ##this is the number of internal knots
degree <- 1
knots <- unname(quantile(backward_time_vector, probs = seq(from = 0, to = 1, length.out = num_knots+2)[-c(1, num_knots+2)]))
#lb <- min(backward_time_vector)
#ub <- max(backward_time_vector)

stan_data <- list(
  N = nrow(final_data),  # Total number of subjects
  K = length(unique(final_data$cluster)),  # Number of clusters
  T = length(grep("time_", names(final_data))),  # Number of time points
  cluster = final_data$cluster,  # Cluster indicator
  treatment = final_data$treatment,
  time_points = time_points,
  x1 = final_data$x1,  # Binary covariate
  x2 = final_data$x2,  # Continuous covariate
  Y = as.matrix(final_data[, grep("time_", names(final_data))]),  # Longitudinal measurements
  MASK = mask,  # Mask for missing values
  survival_time = final_data$observed_time,  # Survival times
  status = final_data$status,  # Censoring indicator
  num_knots = num_knots,
  knots = knots,
  spline_degree = degree
)



# Compile the Stan model once
stan_model_code <- "
functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order == 1) {
      for (i in 1:size(t)) {
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] <= ext_knots[ind + 1]);
      }
    } else {
      if (ext_knots[ind] != ext_knots[ind + order - 1]) {
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind + order - 1] - ext_knots[ind]);
      }
      if (ext_knots[ind + 1] != ext_knots[ind + order]) {
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind + 1], size(t))) /
                 (ext_knots[ind + order] - ext_knots[ind + 1]);
      }
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order - 1) +
                 w2 .* build_b_spline(t, ext_knots, ind + 1, order - 1);
    }
    return b_spline;
  }
}

data {
  int<lower=0> N;  // total number of subjects
  int<lower=0> K;  // number of clusters
  int<lower=0> T;  // number of time points
  int<lower=1,upper=K> cluster[N];  // cluster indicator
  vector[N] x1;  // binary covariate
  vector[N] x2;  // continuous covariate
  vector[N] treatment;
  row_vector[T] time_points;
  matrix[N, T] Y;  // longitudinal measurements
  matrix[N, T] MASK;  // mask for missing values
  vector[N] survival_time;  // observed survival or censoring times
  int<lower=0,upper=1> status[N];  // censoring indicator
  int num_knots;  // number of internal knots for the splines
  vector[num_knots] knots;  // sequence of internal knots
  int spline_degree;  // degree of the splines (order - 1)
}

transformed data {
  matrix[N, T] X1 = rep_matrix(x1, T);
  matrix[N, T] X2 = rep_matrix(x2, T);
  matrix[N, T] TRT = rep_matrix(treatment, T);
  
  int num_basis = num_knots + spline_degree + 1;  // total number of B-splines
  
  int non_missing_count = 0;  // Count of non-missing values
  
  for (i in 1:N) {
    for (t in 1:T) {
      if (MASK[i, t] > 0) {
        non_missing_count += 1;
      }
    }
  }
  
  vector[non_missing_count] Y_non_missing;  // Vector to store non-missing values
  int non_missing_indices[non_missing_count, 2];  // Store (subject, time) pairs of non-missing values
  int index = 1;  // Counter for non-missing values

// Extract non-missing values and their indices
  for (i in 1:N) {
    for (t in 1:T) {
      if (MASK[i, t] > 0) {
        Y_non_missing[index] = Y[i, t];  // Extract non-missing values
        non_missing_indices[index, 1] = i;  // Store the subject index
        non_missing_indices[index, 2] = t;  // Store the time index
        index += 1;
      }
    }
  }

}

parameters {

  real alpha01;
  real alpha02;
  real alpha11;
  real alpha12;

  real b;
  real c;

  real<lower=0> sigma_b;
  real<lower=0> sigma_u;
  real<lower=0> sigma_e;
  real<lower=0> lambda0;
  real<lower=0> gamma;


  row_vector[num_basis] a_backward;  // Raw spline coefficients for backward time
  row_vector[num_basis] a_treatment;  // Raw spline coefficients for treatment effect
  

  vector[N] z_b;
  vector[K] z_u;
  vector<lower=0, upper=1>[N] U;
}

transformed parameters {
  vector[N] lambda;
  vector[N] death_time;
  vector[N] b_i = z_b * sigma_b;
  vector[K] u_i = z_u * sigma_u;

  vector[non_missing_count] backward_times_non_missing;  // Only non-missing backward times
  matrix[num_basis, non_missing_count] B;  
  vector[non_missing_count] spline_contribution_backward;
  vector[non_missing_count] spline_contribution_treatment;

  matrix[N, T] BI = rep_matrix(b_i, T);
  matrix[N, T] UI;
  matrix[N, T] MU;
  
  
  lambda = lambda0 * exp(alpha11 * x1 + alpha12 * x2 + b * b_i + c * u_i[cluster]);

  for (i in 1:N) {
    if (status[i] == 1) {
      death_time[i] = survival_time[i];  // use observed death time for uncensored
    } else {
      death_time[i] = pow(-log(exp(-lambda[i] * pow(survival_time[i], gamma)) - U[i] * exp(-lambda[i] * pow(survival_time[i], gamma))) / lambda[i], 1 / gamma);
    }
  }

  vector[spline_degree + num_knots + 1] ext_knots_temp;
  vector[2 * spline_degree + num_knots + 2] ext_knots;  // extended knots

  // Calculate backward times for non-missing values
  for (ind in 1:non_missing_count) {
    int subj = non_missing_indices[ind, 1];
    int time_idx = non_missing_indices[ind, 2];
    backward_times_non_missing[ind] = death_time[subj] - time_points[time_idx];
  }
  
  // Generate extended knots
  ext_knots_temp = append_row(rep_vector(min(backward_times_non_missing), spline_degree+1), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(max(backward_times_non_missing), spline_degree+1));
  
  // Calculate spline contributions for non-missing values
  for (ind in 1:num_basis) {
    B[ind, :] = to_row_vector(build_b_spline(to_array_1d(backward_times_non_missing), to_array_1d(ext_knots), ind, spline_degree + 1));
  }
  

  spline_contribution_backward = to_vector(a_backward * B);
  spline_contribution_treatment = to_vector(a_treatment * B);
  
  MU = rep_matrix(0, N, T);
  UI = rep_matrix(u_i[cluster], T);
  
  for (ind in 1:non_missing_count) {
    int subj = non_missing_indices[ind, 1];
    int time_idx = non_missing_indices[ind, 2];
    MU[subj, time_idx] = X1[subj, time_idx] * alpha01 + X2[subj, time_idx] * alpha02 + spline_contribution_backward[ind] + TRT[subj, time_idx] * spline_contribution_treatment[ind] + BI[subj, time_idx] + UI[subj, time_idx];
    
  }
  
  

}

model {
  // Priors
  alpha01 ~ normal(0, 5);
  alpha02 ~ normal(0, 5);
  alpha11 ~ normal(0, 5);
  alpha12 ~ normal(0, 5);

  b ~ normal(0, 1);
  c ~ normal(0, 1);

  sigma_b ~ normal(0, 5);
  sigma_u ~ normal(0, 5);
  sigma_e ~ normal(0, 5);
  lambda0 ~ inv_gamma(2, 1);
  gamma ~ gamma(0.5, 0.5);
  
  a_backward ~ normal(0, 15);
  a_treatment ~ normal(0, 15);

  z_b ~ normal(0, 1);
  z_u ~ normal(0, 1);
  U ~ uniform(0, 1);

  // Longitudinal model
   for (i in 1:N) {
    for (t in 1:T) {
      if (MASK[i, t] > 0) {
        Y[i, t] ~ normal(MU[i, t], sigma_e);
      }
    }
   }
  
   // Survival likelihood
  for (i in 1:N) {
    if (status[i] == 1) {
      // Observed death times
      target += weibull_lpdf(survival_time[i] | gamma, pow(lambda[i], -1 / gamma));
    } else {
      // Censored observations
      target += weibull_lccdf(survival_time[i] | gamma, pow(lambda[i], -1 / gamma));
    }
  }
  

}


"

stan_model <- stan_model(model_code = stan_model_code)

init_fn <- function() {
  list(alpha01 = 1, alpha02 = 0.9, alpha11 = 0.2, alpha12 = -0.01,b = 0.03, c = 0.02, sigma_u = 5, sigma_b = 6, sigma_e = 4, lambda0 = 0.05, gamma = 2.2,
       z_b = rnorm(n, 0, 1),
       z_u = rnorm(length(unique(final_data$cluster)), 0, 1),  
       U = runif(n, 0, 1), 
       a_backward = rnorm(length(knots)+degree+1, 0, 1),  
       a_treatment = rnorm(length(knots)+degree+1, 0, 1))
}

# Compile and sample from the Stan model
fit <- sampling(stan_model, data = stan_data, init = init_fn, iter = 2000, warmup = 1000, chains = 2, control = list(adapt_delta = 0.99, max_treedepth = 15), cores=2, refresh=100)


print(fit)  

traceplot(fit, c("alpha01"))


# Extract posterior samples
post_samples <- extract(fit)


backward <- posterior_samples$spline_contribution_backward
treatment<- posterior_samples$spline_contribution_treatment

mean_backward <- apply(backward, 2, mean)
mean_treatment <- apply(treatment, 2, mean)

true_values <- final_data %>%
  tidyr::pivot_longer(cols = starts_with("time_"), names_to = "timepoint", values_to = "measurement") %>%
  mutate(
    timepoint = as.numeric(gsub("time_", "", timepoint)),
    backward_time = survival_time - timepoint,
    True_Backward = alpha00 - (alpha03 / (1 + alpha04 * backward_time)),
    True_Treatment = alpha05 * exp(alpha06 * backward_time + alpha07)
  )

true_values <- true_values %>%
  mutate(
    Estimated_Backward = mean_backward,
    Estimated_Treatment = mean_treatment
  )


ggplot(true_values, aes(x = backward_time)) +
  geom_line(aes(y = True_Backward, color = "True Backward"), linetype = "dashed") +
  geom_line(aes(y = Estimated_Backward, color = "Estimated Backward"), alpha = 0.7) +
  labs(title = "True vs Estimated Spline Contribution (Backward)", x = "Backward Time", y = "Value") +
  scale_color_manual(values = c("True Backward" = "blue", "Estimated Backward" = "red")) +
  theme_minimal()

# Plot the true vs estimated values for treatment spline contribution
ggplot(true_values, aes(x = backward_time)) +
  geom_line(aes(y = True_Treatment, color = "True Treatment"), linetype = "dashed") +
  geom_line(aes(y = Estimated_Treatment, color = "Estimated Treatment"), alpha = 0.7) +
  labs(title = "True vs Estimated Spline Contribution (Treatment)", x = "Backward Time", y = "Value") +
  scale_color_manual(values = c("True Treatment" = "green", "Estimated Treatment" = "orange")) +
  theme_minimal()