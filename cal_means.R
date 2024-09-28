all_means <- list()
path <- "C:/Users/feiyi/OneDrive/Desktop/Katie/Terminal-Decline/Results/"


for (i in 0:49) {
  file_name <- sprintf("%smod0.result.%d.csv", path, i)
  data <- read.csv(file_name)
  all_means[[i + 1]] <- data$mean[1:12]
}


mean_matrix <- do.call(cbind, all_means)
row_means <- rowMeans(mean_matrix, na.rm = TRUE)

parameter_names <- read.csv(sprintf("%smod0.result.0.csv", path))$X[1:12]
output_df <- data.frame(Parameter = parameter_names, Mean = row_means)
write.csv(output_df, file = sprintf("%spar_means_output.csv", path), row.names = FALSE)