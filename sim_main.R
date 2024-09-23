#!/usr/bin/env Rscript
library(MASS)
library(survival)
library(tidyr)
library(dplyr)
library(rstan)

mask <- as.matrix(read.csv(list.files(pattern="mask.")))
time_points <- c(0, 3, 6, 9, 12)
final_data <- as.data.frame(read.csv(list.files(pattern="sim.data.")))

  stan_data <- list(
    N = nrow(final_data),  # Total number of subjects
    K = length(unique(final_data$cluster)),  # Number of clusters
    T = length(grep("time_", names(final_data))),  # Number of time points
    cluster = final_data$cluster,  # Cluster indicator
    treatment = final_data$treatment,
    x1 = final_data$x1,  # Binary covariate
    x2 = final_data$x2,  # Continuous covariate
    y = as.matrix(final_data[, grep("time_", names(final_data))]),  # Longitudinal measurements
    mask = mask,  # Mask for missing values
    survival_time = final_data$time,  # Survival times
    status = final_data$status,  # Censoring indicator
    time_points = time_points
  )
  
  
  
  # Compile the Stan model once
  stan_model_code <- "
data {
  int<lower=0> N;  // total number of subjects
  int<lower=0> K;  // number of clusters
  int<lower=0> T;  // number of time points
  int<lower=1,upper=K> cluster[N];  // cluster indicator
  vector[N] x1;  // binary covariate
  vector[N] x2;  // continuous covariate
  vector[N] treatment;
  matrix[N, T] y;  // longitudinal measurements
  matrix[N, T] mask;  // mask for missing values
  real<lower=0> time_points[T];
  real<lower=0> survival_time[N];  // observed survival or censoring times
  int<lower=0,upper=1> status[N];  // censoring indicator
}

parameters {
  real alpha00;
  real alpha01;
  real alpha02;
  real alpha03;
  real alpha04;
  real alpha11;
  real alpha12;
  real alpha13;
  real b;
  real c;
  real<lower=0> sigma_b;
  real<lower=0> sigma_u;
  real<lower=0> sigma_e;
  real<lower=0> lambda0;
  real<lower=0> gamma;
  
  vector[N] z_b;  
  vector[K] z_u;
  vector<lower=0,upper=1>[N] U;
}

transformed parameters {
  
  vector[N] b_i = z_b * sigma_b; 
  vector[K] u_i = z_u * sigma_u;
  vector[N] death_time;
  vector[N] F_C;
  vector[N] U_adjusted;
  vector[N] eta;

  for (i in 1:N) {
    eta[i] = alpha11 * x1[i] + alpha12 * x2[i] + alpha13 * treatment[i] + c * u_i[cluster[i]] + b * b_i[i];
  }
   
  
  for (i in 1:N) {
    if (status[i] == 1) {
      death_time[i] = survival_time[i];  // use observed death time for uncensored
    } else {
      F_C[i] = 1 - exp(-(lambda0 * exp(eta[i])) * pow(survival_time[i], gamma));
      U_adjusted[i] = F_C[i] + U[i] * (1 - F_C[i]);
      death_time[i] = pow(-log(1 - U_adjusted[i]) / (lambda0 * exp(eta[i])), 1 / gamma);
    }
  }
}

model {
  // Priors
  alpha00 ~ normal(0, 10);
  alpha01 ~ normal(0, 10);
  alpha02 ~ normal(0, 10);
  alpha03 ~ normal(0, 10);
  alpha04 ~ normal(0, 10);
  alpha11 ~ normal(0, 10);
  alpha12 ~ normal(0, 10);
  alpha13 ~ normal(0, 10);
  b ~ normal(0, 10);
  c ~ normal(0, 10);

  sigma_b ~ normal(0, 10);
  sigma_u ~ normal(0, 10);
  sigma_e ~ normal(0, 10);
  lambda0 ~ inv_gamma(0.01, 0.01);
  gamma ~ gamma(0.01, 0.01);

  z_b ~ normal(0, 1);
  z_u ~ normal(0, 1);
  
  

  // Longitudinal model
  for (i in 1:N) {
    for (t in 1:T) {
      if (mask[i, t]) {
        real back_t = death_time[i] - time_points[t]
        y[i, t] ~ normal(alpha00 + x1[i] * alpha01 + x2[i] * alpha02 + back_t * alpha03 + treatment[i] * alpha04 + b_i[i] + u_i[cluster[i]], sigma_e);
      }
    }
  }
  
    
  // Survival likelihood
 for (i in 1:N) {
    if (status[i] == 1) {
      // Observed death time: use weibull exponential likelihood
      target += log(lambda0) + log(gamma) + (gamma - 1) * log(death_time[i]) + eta[i] - lambda0 * exp(eta[i]) * death_time[i]^gamma;
    } else {
      // Censored event: Weibull survival function
      target += -lambda0 * exp(eta[i]) * survival_time[i]^gamma;
    }
  }
}
"
stan_model <- stan_model(model_code = stan_model_code)

init_fn <- function() {
  list(alpha00 = 10, alpha01 = 5, alpha02 = 0.3, alpha03 = 2, alpha04 = 10, alpha11 = -2, alpha12 = 0.02, alpha13 = -4, b = 0.5, c = 0.4, lambda0 = 0.03, gamma = 1.5, sigma_b = 1, sigma_u = 3, sigma_e = 3)
}

# Compile and sample from the Stan model
fit <- sampling(stan_model, data = stan_data, init = init_fn, iter = 2000, warmup = 1000, chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 15), cores=4)

result <- summary(fit)
fit_df <- as.data.frame(result$summary)
text <- list.files(pattern="sim.data.")
num <- unlist(lapply(strsplit(text,'.',fixed=TRUE),function(x) x[[3]]))
write.csv(fit_df, paste0("mod.result.",num,".csv"))


pdf(file = paste0("mod.traceplot.",num,".pdf"),   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
traceplot(fit, c("alpha00","alpha01","alpha02","alpha03","alpha04","alpha11","alpha12","alpha13","b","c","lambda0","gamma","sigma_b","sigma_u","sigma_e"))
dev.off()


