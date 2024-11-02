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
quan_len<-seq(from=0,to=1,length.out = num_knots+2)
knots<-as.numeric(quantile(backward_time_vector,quan_len))

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
  knots = knots
)


stan_model_code <- "
data {
  int<lower=0> N;  
  int<lower=0> K;  
  int<lower=0> T;  
  int<lower=1,upper=K> cluster[N];  
  vector[N] x1;  
  vector[N] x2;  
  vector[N] treatment;
  row_vector[T] time_points;
  matrix[N, T] Y;  
  matrix[N, T] MASK;  
  vector[N] survival_time;  
  int<lower=0,upper=1> status[N];  
  int num_knots;
  vector[num_knots+2] knots;  
}

transformed data {
  // Count the total number of non-missing values
  int non_missing_count = 0;
  for (i in 1:N) {
    for (t in 1:T) {
      if (MASK[i, t] > 0) {
        non_missing_count += 1;
      }
    }
  }
  
  // Store indices of non-missing values
  int non_missing_indices[non_missing_count, 2];
  int index = 1;
  for (i in 1:N) {
    for (t in 1:T) {
      if (MASK[i, t] > 0) {
        non_missing_indices[index, 1] = i;
        non_missing_indices[index, 2] = t;
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
  vector[N] z_b;
  vector[K] z_u;
  vector<lower=0, upper=1>[N] U;
  vector[num_knots+1] beta_piecewise_0; 
  vector[num_knots+1] beta_piecewise_1; 
}


transformed parameters {
  vector[N] lambda;
  vector[N] death_time;
  vector[N] b_i = z_b * sigma_b;
  vector[K] u_i = z_u * sigma_u;
  matrix[non_missing_count, num_knots+1] piecewise_matrix;
  vector[non_missing_count] MU_non_missing;

  lambda = lambda0 * exp(alpha11 * x1 + alpha12 * x2 + b * b_i + c * u_i[cluster]);

  for (i in 1:N) {
    if (status[i] == 1) {
      death_time[i] = survival_time[i];  // Observed death time
    } else {
      death_time[i] = pow(-log(exp(-lambda[i] * pow(survival_time[i], gamma)) - U[i] * exp(-lambda[i] * pow(survival_time[i], gamma))) / lambda[i], 1 / gamma);
    }
  }

  // Determine maximum death time across all subjects
  real max_death_time = max(death_time);

  // Extend knots to include the maximum death time
  vector[num_knots + 2] extended_knots;
  for (k in 2:num_knots+1) {
    extended_knots[k] = knots[k];
  }
  extended_knots[1] = 0;
  extended_knots[num_knots + 2] = max_death_time;

  // Calculate piecewise matrix based on death_time for non-missing values only
  for (ind in 1:non_missing_count) {
    int i = non_missing_indices[ind, 1];
    int t = non_missing_indices[ind, 2];
    real backward_time = death_time[i] - time_points[t];
    for (k in 1:num_knots+1) {
      if (backward_time > extended_knots[k] && backward_time <= extended_knots[k + 1]) {
        piecewise_matrix[ind, k] = backward_time - extended_knots[k];
      } else if (backward_time > extended_knots[k + 1]) {
        piecewise_matrix[ind, k] = extended_knots[k + 1] - extended_knots[k];
      } else {
        piecewise_matrix[ind, k] = 0;
      }
    }
  }

  // Calculate MU for non-missing values
  for (ind in 1:non_missing_count) {
    int i = non_missing_indices[ind, 1];
    MU_non_missing[ind] = alpha01 * x1[i] + alpha02 * x2[i] + dot_product(piecewise_matrix[ind, ], beta_piecewise_0) +
                          treatment[i] * dot_product(piecewise_matrix[ind, ], beta_piecewise_1) +
                          b_i[i] + u_i[cluster[i]];
  }
}

model {
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
  z_b ~ normal(0, 1);
  z_u ~ normal(0, 1);
  U ~ uniform(0, 1);
  beta_piecewise_0 ~ normal(0, 5); 
  beta_piecewise_1 ~ normal(0, 5); 

  // Longitudinal likelihood for non-missing values
  for (ind in 1:non_missing_count) {
    Y[non_missing_indices[ind, 1], non_missing_indices[ind, 2]] ~ normal(MU_non_missing[ind], sigma_e);
  }

  // Survival likelihood
  for (i in 1:N) {
    if (status[i] == 1) {
      target += weibull_lpdf(survival_time[i] | gamma, pow(lambda[i], -1 / gamma));
    } else {
      target += weibull_lccdf(survival_time[i] | gamma, pow(lambda[i], -1 / gamma));
    }
  }
}
"

stan_model <- stan_model(model_code = stan_model_code)

init_fn <- function() {
  list(alpha01 = 5, alpha02 = 0.3, alpha11 = -1, alpha12 = 0.05,b = 0.03, c = 0.05, sigma_u = 5, sigma_b = 8, sigma_e = 3, lambda0 = 0.03, gamma = 1.8,
       z_b = rnorm(n, 0, 1),
       z_u = rnorm(length(unique(final_data$cluster)), 0, 1),  
       U = runif(n, 0, 1), 
       beta_piecewise_0 = rnorm(3, 0, 1),  
       beta_piecewise_1 = rnorm(3, 0, 1))
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
traceplot(fit, c("alpha01","alpha02","alpha11","alpha12","b","c","sigma_u","sigma_b", "sigma_e"))
dev.off()


