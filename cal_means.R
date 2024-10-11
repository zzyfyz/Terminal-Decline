all_means <- list()
all_sds <- list()
coverage_counts <- rep(0, 14)

path <- "C:/Users/feiyi/OneDrive/Desktop/Katie/Terminal-Decline/Results/"

# True values of the parameters
true_values <- c(
  alpha00 = 10,
  alpha01 = 5,
  alpha02 = 0.3,
  alpha03 = 2,
  alpha04 = 10,
  alpha11 = -1,
  alpha12 = 0.05,
  b = 0.03,
  c = 0.05,
  sigma_b = 8,
  sigma_u = 5, 
  sigma_e = 3,
  lambda0 = 0.04,
  gamma = 1.2
)

for (i in 0:99) {
  file_name <- sprintf("%smod.result.%d.csv", path, i)
  data <- read.csv(file_name)
  
  # Store the means and standard deviations
  all_means[[i + 1]] <- data$mean[1:14]
  all_sds[[i + 1]] <- data$sd[1:14]
  
  # Calculate coverage probability for each parameter
  for (j in 1:14) {
    if (true_values[j] >= data$X2.5.[j] && true_values[j] <= data$X97.5.[j]) {
      coverage_counts[j] <- coverage_counts[j] + 1
    }
  }
}

mean_matrix <- do.call(cbind, all_means)
sd_matrix <- do.call(cbind, all_sds)

row_means <- rowMeans(mean_matrix, na.rm = TRUE)
row_sd_means <- rowMeans(sd_matrix, na.rm = TRUE)

# Calculate coverage probabilities
coverage_probabilities <- coverage_counts / 100

# Read parameter names from the first file
parameter_names <- read.csv(sprintf("%smod.result.0.csv", path))$X[1:14]
output_df <- data.frame(
  Parameter = parameter_names,
  True_Value = true_values,
  Mean = row_means,
  Mean_SD = row_sd_means,
  Coverage_Probability = coverage_probabilities
)
write.csv(output_df, file = sprintf("%spar_means_output.csv", path), row.names = FALSE)