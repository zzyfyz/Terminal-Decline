library(readr)
library(dplyr)
library(stringr)
library(purrr)


# Define true values for the parameters
true_values <- c(alpha1 = 0.7, alpha2 = -0.03, phi0 = -4.8, phi1 = -0.02, phi2 = 0.01, v1 = 0.3, v2 = -0.2, 
                 sigma_b = 1, sigma_u = 0.5, 
                 sigma_e = 1, gamma = 1.4)

# Set the directory path
file_path <- "C:/Yizhou/Results/C10_D40_512"

files <- list.files(path = file_path, pattern = "^mod\\.result\\.(\\d+)\\.csv$", full.names = TRUE)
idx   <- as.integer(str_match(basename(files), "^mod\\.result\\.(\\d+)\\.csv$")[,2])
o     <- order(idx); files <- files[o]

read_one <- function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  
  # rename first column to Parameter (handles ...1 or already-named exports)
  if ("Parameter" %in% names(df)) {
    # nothing
  } else if ("...1" %in% names(df)) {
    df <- df %>% rename(Parameter = `...1`)
  } else {
    stop(sprintf("Can't find parameter-name column in %s", f))
  }
  
  df %>%
    mutate(Parameter = trimws(Parameter)) %>%
    filter(Parameter %in% names(true_values)) %>%
    transmute(
      file  = basename(f),
      Parameter,
      est   = mean,
      lower = `2.5%`,
      upper = `97.5%`,
      true  = unname(true_values[Parameter]),
      bias  = est - true,
      mse   = (est - true)^2,
      coverage = as.integer(true >= lower & true <= upper)
    )
}

all_results <- map_dfr(files, read_one)

results_summary <- all_results %>%
  group_by(Parameter) %>%
  summarise(
    True_Value = first(true),
    Bias       = mean(bias, na.rm = TRUE),
    MSE        = mean(mse,  na.rm = TRUE),
    Coverage   = mean(coverage, na.rm = TRUE),
    .groups = "drop"
  )

# Save the summary to a CSV file
write.csv(results_summary, file = file.path(file_path, "simulation_summary.csv"), row.names = FALSE)

# Print the summary
print(results_summary)
