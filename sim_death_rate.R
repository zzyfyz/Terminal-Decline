# Define how many datasets and where they are located
num_datasets <- 300
data_path <- "C:/Yizhou/Term/Terminal-Decline/Spline/"  # e.g., "C:/Users/YourName/Desktop/"

# Initialize a vector to store the death rate for each dataset
death_rates <- numeric(num_datasets)

# Loop over the 300 datasets
for (i in 0:(num_datasets - 1)) {
  # Construct the filename (assuming files are named sim.data.0.csv, sim.data.1.csv, etc.)
  file_name <- paste0(data_path, "sim.data.", i, ".csv")
  
  # Read the CSV file
  dat <- read.csv(file_name)
  
  # Calculate the death rate as the proportion of 'status' == 1
  # (assuming 'status' is 1 for death, 0 for censored/surviving)
  death_rate <- mean(dat$status == 1, na.rm = TRUE)
  
  # Store the result
  death_rates[i + 1] <- death_rate
}

# Examine the results
summary(death_rates)
