dirg <- "C:/Yizhou/Term/Terminal-Decline/Piecewise"  

files <- list.files(
  path = dirg,
  pattern = "^sim\\.data\\.(\\d+)\\.csv$",
  full.names = TRUE
)

# sort by numeric suffix
idx <- as.integer(sub(".+sim\\.data\\.(\\d+)\\.csv$", "\\1", files))
o <- order(idx); files <- files[o]; idx <- idx[o]

death_rates <- sapply(files, function(f) {
  df <- read.csv(f)
  mean(df$status == 1, na.rm = TRUE)
})

avg_death_rate <- mean(death_rates)
sd_death_rate  <- sd(death_rates)

pooled_counts <- sapply(files, function(f) {
  df <- read.csv(f)
  c(deaths = sum(df$status == 1, na.rm = TRUE),
    n      = sum(!is.na(df$status)))
})
pooled_death_rate <- sum(pooled_counts["deaths", ]) / sum(pooled_counts["n", ])

# save next to the data
out <- data.frame(file = basename(files), index = idx, death_rate = as.numeric(death_rates))
write.csv(out, file.path(dirg, "death_rates_by_dataset.csv"), row.names = FALSE)
