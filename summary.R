all_means <- list()
coverage_counts <- rep(0, 11)

path <- "C:/Yizhou/Results/Piecewise_412"

# True values (make sure order matches your result files)
true_values <- c(
  alpha1 = 0.7,
  alpha2 = -0.03,
  phi0   = -4.7,
  phi1   = -0.02,
  phi2   = 0.01,
  v1     = 0.3,
  v2     = -0.2,
  sigma_b= 1,
  sigma_u= 0.5, 
  sigma_e= 1,
  gamma  = 1.4
)

# Get all matching files
file_list <- list.files(path, pattern = "mod\\.result\\..*\\.csv$", full.names = TRUE)
num_files <- length(file_list)
if (num_files == 0) stop("No result files found.")

# Pre-allocate matrices for speed/clarity
K <- 11  # number of parameters summarized per file
mean_matrix <- matrix(NA_real_, nrow = K, ncol = num_files)
low_matrix  <- matrix(NA_real_, nrow = K, ncol = num_files)
high_matrix <- matrix(NA_real_, nrow = K, ncol = num_files)

# Read first file to determine parameter names robustly
first_df <- read.csv(file_list[1], check.names = FALSE)
parameter_names <- first_df[[1]][1:K]

# Fill matrices & coverage counts
for (i in seq_along(file_list)) {
  df <- read.csv(file_list[i], check.names = FALSE)
  mean_matrix[, i] <- as.numeric(df$mean[1:K])
  # Support either "2.5%" or "X2.5." column names
  low_col  <- if ("2.5%"  %in% names(df)) "2.5%"  else if ("X2.5."  %in% names(df)) "X2.5."  else stop("No 2.5% column.")
  high_col <- if ("97.5%" %in% names(df)) "97.5%" else if ("X97.5." %in% names(df)) "X97.5." else stop("No 97.5% column.")
  low_matrix[, i]  <- as.numeric(df[[low_col]][1:K])
  high_matrix[, i] <- as.numeric(df[[high_col]][1:K])
  
  # Coverage
  coverage_counts <- coverage_counts +
    as.integer((true_values[1:K] >= low_matrix[, i]) & (true_values[1:K] <= high_matrix[, i]))
}

# Row summaries
row_means <- rowMeans(mean_matrix, na.rm = TRUE)
sq_err    <- (mean_matrix - matrix(true_values[1:K], nrow = K, ncol = num_files))^2
mse_values <- rowMeans(sq_err, na.rm = TRUE)

# MCSEs
n <- num_files
row_sds <- apply(mean_matrix, 1, sd, na.rm = TRUE)
bias_values <- row_means - true_values[1:K]
bias_mcse   <- row_sds / sqrt(n)

mse_sd   <- apply(sq_err, 1, sd, na.rm = TRUE)
mse_mcse <- mse_sd / sqrt(n)

coverage_probabilities <- coverage_counts / n
coverage_mcse <- sqrt(coverage_probabilities * (1 - coverage_probabilities) / n)

# Assemble output
output_df <- data.frame(
  Parameter = parameter_names,
  True_Value = as.numeric(true_values[1:K]),
  Mean = row_means,
  Bias = bias_values,
  Bias_MCSE = bias_mcse,
  MSE = mse_values,
  MSE_MCSE = mse_mcse,
  Coverage_Probability = coverage_probabilities,
  Coverage_MCSE = coverage_mcse,
  stringsAsFactors = FALSE
)

# Write output
out_file <- file.path(path, "par_means_output.csv")
write.csv(output_df, file = out_file, row.names = FALSE)