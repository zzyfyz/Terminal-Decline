#!/usr/bin/env Rscript
library(MASS)
library(survival)
library(tidyr)
library(dplyr)
library(rstan)

mask <- as.matrix(read.csv(list.files(pattern="mask.")))
final_data <- as.data.frame(read.csv(list.files(pattern="sim.data.")))
time_points <- seq(0, 12, by = 3)

stan_data <- list(
  N = nrow(final_data),  # Total number of subjects
  K = length(unique(final_data$cluster)),  # Number of clusters
  T = length(grep("time_", names(final_data))),  # Number of time points
  cluster = final_data$cluster,  # Cluster indicator
  treatment = final_data$treatment,
  time_points = time_points,
  x1 = final_data$x1,  # Binary covariate
  x2 = final_data$x2,  # Continuous covariate
  y = as.matrix(final_data[, grep("time_", names(final_data))]),  # Longitudinal measurements
  mask = mask,  # Mask for missing values
  survival_time = final_data$observed_time,  # Survival times
  status = final_data$status  # Censoring indicator
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
  vector[T] time_points;
  matrix[N, T] y;  // longitudinal measurements
  matrix[N, T] mask;  // mask for missing values
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
  
  vector[N] lambda;
  vector[N] death_time;
  vector[N]  b_i = z_b * sigma_b;
  vector[K]  u_i = z_u * sigma_u;
  
  for (i in 1:N) {
    lambda[i] = lambda0 * exp(alpha11 * x1[i] + alpha12 * x2[i] + c * u_i[cluster[i]] + b * b_i[i]);
  }
   
  
  for (i in 1:N) {
    if (status[i] == 1) {
      death_time[i] = survival_time[i];  // use observed death time for uncensored
    } else {
      death_time[i] = pow(-log(exp(-lambda[i] * pow(survival_time[i], gamma))-U[i] *exp(-lambda[i] * pow(survival_time[i], gamma))) / lambda[i], 1 / gamma);
    }
  }
}

model {

  // Priors
  alpha00 ~ normal(0, 5);
  alpha01 ~ normal(0, 5);
  alpha02 ~ normal(0, 5);
  alpha03 ~ normal(0, 5);
  alpha04 ~ normal(0, 5);
  alpha11 ~ normal(0, 5);
  alpha12 ~ normal(0, 5);

  b ~ normal(0, 1);
  c ~ normal(0, 1);

  sigma_b ~ normal(0, 5);
  sigma_u ~ normal(0, 5);
  sigma_e ~ normal(0, 5);
  lambda0 ~ inv_gamma(0.01, 0.01);
  gamma ~ gamma(0.01, 0.01);

  z_b ~ normal(0, 1);
  z_u ~ normal(0, 1);
  U ~ uniform(0,1);
  

  // Longitudinal model
  for (i in 1:N) {
    for (t in 1:T) {
      real backward_time = death_time[i] - time_points[t];
      if (mask[i, t] > 0) {
        y[i, t] ~ normal(alpha00 + x1[i] * alpha01 + x2[i] * alpha02 + backward_time * alpha03 + treatment[i] * alpha04 + b_i[i] + u_i[cluster[i]], sigma_e);
      }
    }
  }
  
    
  // Survival likelihood
 for (i in 1:N) {
    if (status[i] == 1) {
      // Observed death time: use weibull exponential likelihood
      target += log(lambda[i]) + log(gamma) + (gamma - 1) * log(death_time[i])  - lambda[i] * death_time[i]^gamma;
    } else {
      // Censored event: Weibull survival function
      target += -lambda[i] * survival_time[i]^gamma;
    }
  }
}
"
stan_model <- stan_model(model_code = stan_model_code)

init_fn <- function() {
  list(alpha00 = 10, alpha01 = 5, alpha02 = 0.3, alpha03 = 2, alpha04 = 10, alpha11 = 1, alpha12 = 0.05, b = 0.03, c = 1, lambda0 = 0.05, gamma = 1.2, sigma_u = 1, sigma_b = 8, sigma_e = 3)
}

# Compile and sample from the Stan model
fit <- sampling(stan_model, data = stan_data, init = init_fn, iter = 4000, warmup = 2000, chains = 2, control = list(adapt_delta = 0.99, max_treedepth = 15), cores=2)

result <- summary(fit)
fit_df <- as.data.frame(result$summary)
text <- list.files(pattern="sim.data.")
num <- unlist(lapply(strsplit(text,'.',fixed=TRUE),function(x) x[[3]]))
write.csv(fit_df, paste0("mod.result.",num,".csv"))


pdf(file = paste0("mod.traceplot.",num,".pdf"),   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 8) # The height of the plot in inches
traceplot(fit, c("alpha00","alpha01","alpha02","alpha03","alpha04","alpha11","alpha12","b","c","lambda0","gamma","sigma_b","sigma_u","sigma_e"))
dev.off()


