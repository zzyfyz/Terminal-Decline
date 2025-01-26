library(readr)
library(dplyr)
library(stringr)

# Define true values for the parameters
true_values <- c(alpha01 = 1, alpha02 = 1.6, alpha11 = 0.2, alpha12 = -0.01, 
                 b = 0.03, c = 0.02, sigma_b = 6, sigma_u = 5, 
                 sigma_e = 4, lambda0 = 0.05, gamma = 2.2)

# Set the directory path
file_path <- "C:/Yizhou/Results"

# Get all simulation result files
file_list <- list.files(path = file_path, pattern = "mod\\.result\\..*\\.csv", full.names = TRUE)

# Initialize a data frame to store the summary
results_summary <- data.frame(Parameter = names(true_values),
                              True_Value = true_values,
                              Bias = NA,
                              MSE = NA,
                              Coverage = NA)

# Loop over each file and summarize results
all_results <- lapply(file_list, function(file) {
  # Read the file
  data <- read_csv(file, col_names = TRUE)
  
  # Assign proper column names based on the structure
  colnames(data) <- c("Parameter", "mean", "se_mean", "sd", "2.50%", "25%", "50%", "75%", "97.50%", "n_eff", "Rhat")
  
  # Calculate metrics for each parameter
  metrics <- lapply(names(true_values), function(param) {
    true_value <- true_values[param]
    
    # Filter data for the current parameter
    param_data <- data %>% filter(Parameter  == param)
    
    # Extract relevant columns
    est <- param_data$mean
    lower <- param_data$`2.50%`
    upper <- param_data$`97.50%`
    
    # Calculate bias, MSE, and 95% coverage
    bias <- mean(est - true_value)
    mse <- mean((est - true_value)^2)
    coverage <- mean(true_value >= lower & true_value <= upper)
    
    return(c(bias = bias, mse = mse, coverage = coverage))
  })
  
  # Combine results for all parameters in this file
  do.call(rbind, metrics)
})

# Combine results across all files
all_results_df <- do.call(rbind, all_results)

# Summarize results across all simulations
for (i in seq_along(names(true_values))) {
  param <- names(true_values)[i]
  param_results <- all_results_df[param, ]
  
  results_summary[i, "Bias"] <- mean(param_results[, "bias"])
  results_summary[i, "MSE"] <- mean(param_results[, "mse"])
  results_summary[i, "Coverage"] <- mean(param_results[, "coverage"])
}

# Save the summary to a CSV file
write.csv(results_summary, file = file.path(file_path, "simulation_summary.csv"), row.names = FALSE)

# Print the summary
print(results_summary)
