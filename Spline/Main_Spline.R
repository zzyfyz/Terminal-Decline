#!/usr/bin/env Rscript
library(MASS)
library(survival)
library(tidyr)
library(dplyr)
library(rstan)

mask <- as.matrix(read.csv(list.files(pattern="mask.")))
final_data <- as.data.frame(read.csv(list.files(pattern="sim.data.")))
time_points <- seq(0, 12, by = 3)
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

num_knots <- 2
knots <- unname(quantile(backward_time_vector, probs = seq(from = 0.01, to = 0.99, length.out = num_knots)))

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
  spline_degree = 3
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
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind + 1]);
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
  vector[num_knots] knots;  // sequence of knots
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


  real<lower=0> sigma_e;


  row_vector[num_basis] a_raw_backward;  // Raw spline coefficients for backward time
  row_vector[num_basis] a_raw_treatment;  // Raw spline coefficients for treatment effect
  real<lower=0> tau_backward;  // Step size for the random walk (backward time)
  real<lower=0> tau_treatment;  // Step size for the random walk (treatment effect)


  vector<lower=0, upper=1>[N] U;
}

transformed parameters {
  vector[N] lambda;
  vector[N] death_time;

  vector[non_missing_count] backward_times_non_missing;  // Only non-missing backward times
  matrix[num_basis, non_missing_count] B;  
  vector[non_missing_count] spline_contribution_backward;
  vector[non_missing_count] spline_contribution_treatment;

  row_vector[num_basis] a_backward; 
  row_vector[num_basis] a_treatment; 


  matrix[N, T] MU;


  a_backward = a_raw_backward*tau_backward;
  a_treatment = a_raw_treatment*tau_treatment;
  
  
  
  lambda = 0.03 * exp(alpha11 * x1 + alpha12 * x2 );

  for (i in 1:N) {
    if (status[i] == 1) {
      death_time[i] = survival_time[i];  // use observed death time for uncensored
    } else {
      death_time[i] = pow(-log(exp(-lambda[i] * pow(survival_time[i], 1.8)) - U[i] * exp(-lambda[i] * pow(survival_time[i], 1.8))) / lambda[i], 1 / 1.8);
    }
  }

  vector[spline_degree + num_knots + 1] ext_knots_temp;
  vector[2 * spline_degree + num_knots + 2] ext_knots;  // extended knots

  // Generate extended knots
  ext_knots_temp = append_row(rep_vector(min(death_time), spline_degree+1), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(max(death_time), spline_degree+1));

  // Calculate backward times for non-missing values
  for (ind in 1:non_missing_count) {
    int subj = non_missing_indices[ind, 1];
    int time_idx = non_missing_indices[ind, 2];
    backward_times_non_missing[ind] = death_time[subj] - time_points[time_idx];
  }
  
  // Calculate spline contributions for non-missing values
  for (ind in 1:num_basis) {
    B[ind, :] = to_row_vector(build_b_spline(to_array_1d(backward_times_non_missing), to_array_1d(ext_knots), ind, spline_degree + 1));
  }
  

  spline_contribution_backward = to_vector(a_backward * B);
  spline_contribution_treatment = to_vector(a_treatment * B);
  
  MU = rep_matrix(0, N, T);
  
  for (ind in 1:non_missing_count) {
    int subj = non_missing_indices[ind, 1];
    int time_idx = non_missing_indices[ind, 2];
    MU[subj, time_idx] = X1[subj, time_idx] * alpha01 + X2[subj, time_idx] * alpha02 + spline_contribution_backward[ind] + TRT[subj, time_idx] * spline_contribution_treatment[ind];
    
  }
  
  

}

model {
  // Priors
  alpha01 ~ normal(0, 5);
  alpha02 ~ normal(0, 5);
  alpha11 ~ normal(0, 5);
  alpha12 ~ normal(0, 5);



  sigma_e ~ normal(0, 5);
  
  tau_backward ~ normal(0, 1);
  tau_treatment ~ normal(0, 1);
  a_raw_backward ~ normal(0, 1);
  a_raw_treatment ~ normal(0, 1);


  U ~ uniform(0, 1);

  // Longitudinal model
   for (i in 1:N) {
    for (t in 1:T) {
      if (MASK[i, t] > 0) {
        Y[i, t] ~ normal(MU[i, t], sigma_e);
      }
    }
   }
  

}


"

stan_model <- stan_model(model_code = stan_model_code)

init_fn <- function() {
  list(alpha01 = 5, alpha02 = 0.3, alpha11 = -1, alpha12 = 0.05,sigma_e = 3, 
       tau_backward = 1,
       tau_treatment = 1)
}

# Compile and sample from the Stan model
fit <- sampling(stan_model, data = stan_data, init = init_fn, iter = 2000, warmup = 1000, chains = 2, control = list(adapt_delta = 0.99, max_treedepth = 15), cores=2)

result <- summary(fit)
fit_df <- as.data.frame(result$summary)
text <- list.files(pattern="sim.data.")
num <- unlist(lapply(strsplit(text,'.',fixed=TRUE),function(x) x[[3]]))
write.csv(fit_df, paste0("mod.result.",num,".csv"))


pdf(file = paste0("mod.traceplot.",num,".pdf"),   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 8) # The height of the plot in inches
traceplot(fit, c("alpha01","alpha02","alpha11","alpha12","sigma_e"))
dev.off()


