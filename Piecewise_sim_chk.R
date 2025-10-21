library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# --- Inputs ---
dirg    <- "C:/Yizhou/Sim_data/Piecewise_124"
breaks  <- c(0, 4, 8, 12)       # administrative cap at 12 in your sim
cutlabs <- paste0("(", head(breaks, -1), ",", tail(breaks, -1), "]")

# --- Read & pool all simulated datasets ---
files <- list.files(path = dirg, pattern = "^sim\\.data\\.(\\d+)\\.csv$", full.names = TRUE)
idx   <- as.integer(sub(".+sim\\.data\\.(\\d+)\\.csv$", "\\1", files))
o     <- order(idx); files <- files[o]; idx <- idx[o]

dat_list <- map2(files, idx, ~{
  df <- read.csv(.x)
  df$sim_id <- .y
  df
})

dat <- bind_rows(dat_list)

# Keep only the columns we need
dat <- dat %>%
  transmute(sim_id, time = pmin(observed_time, max(breaks)), status = as.integer(status))

# =============================
# 1) Nelson–Aalen cumulative hazard (pooled)
# =============================
sf <- survfit(Surv(time, status) ~ 1, data = dat)
na_curve <- tibble(
  time = sf$time,
  H_na = -log(sf$surv)  # Nelson–Aalen approx via -log(KM)
)

# =============================
# 2) Piecewise constant hazard via Poisson regression (pooled)
# =============================
# Split follow-up by the pre-specified breaks and estimate interval-specific hazards
split_dat <- survSplit(Surv(time, status) ~ 1, data = dat,
                       cut = breaks[-1], # internal cutpoints
                       episode = "interval",
                       start = "tstart", end = "tstop")

# Compute time at risk and event indicator in each interval row
split_dat <- split_dat %>%
  mutate(
    interval = factor(interval, labels = cutlabs),
    y = status,
    offset_log_pt = log(pmax(tstop - tstart, .Machine$double.eps))
  )

# Poisson model with interval indicators -> hazard per interval
fit_poi <- glm(y ~ interval - 1 + offset(offset_log_pt),
               data = split_dat, family = poisson())

haz_tbl <- tibble(
  interval = levels(split_dat$interval),
  h_hat = exp(coef(fit_poi))
)

# Build a piecewise-linear cumulative hazard from the interval hazards
seg_tbl <- tibble(
  t_left  = head(breaks, -1),
  t_right = tail(breaks, -1)
) %>%
  mutate(interval = factor(paste0("(", t_left, ",", t_right, "]"), levels = cutlabs)) %>%
  left_join(haz_tbl, by = "interval") %>%
  mutate(width = t_right - t_left,
         H_incr = h_hat * width)

seg_tbl <- seg_tbl %>%
  mutate(H_cum_start = c(0, head(cumsum(H_incr), -1)),
         H_cum_end   = cumsum(H_incr))

# Create a stepwise curve for plotting
pw_curve <- bind_rows(
  seg_tbl %>% transmute(time = t_left,  H_pw = H_cum_start),
  seg_tbl %>% transmute(time = t_right, H_pw = H_cum_end)
) %>% arrange(time) %>% distinct(time, .keep_all = TRUE)

# =============================
# 3) Plots
# =============================

# (A) Cumulative hazard: Nelson–Aalen vs piecewise-constant estimate
ggplot() +
  geom_step(data = na_curve, aes(x = time, y = H_na), linewidth = 0.8) +
  geom_line(data = pw_curve, aes(x = time, y = H_pw), linewidth = 0.9, linetype = 2) +
  geom_vline(xintercept = breaks, alpha = 0.4) +
  labs(title = "Cumulative Hazard: Nelson–Aalen vs Piecewise-Constant Estimate (Pooled)",
       x = "Time (months)", y = "Cumulative hazard") +
  theme_minimal(base_size = 14)

# (B) KM survival with cutpoint guides
sf_km <- survfit(Surv(time, status) ~ 1, data = dat)
km_df <- tibble(time = sf_km$time, S = sf_km$surv)

ggplot(km_df, aes(time, S)) +
  geom_step(linewidth = 0.8) +
  geom_vline(xintercept = breaks, alpha = 0.4) +
  labs(title = "Kaplan–Meier Survival with Piecewise Cutpoints (Pooled)",
       x = "Time (months)", y = "S(t)") +
  theme_minimal(base_size = 14)

# =============================
# 4) Optional: per-simulation check (just switch sim_id)
# =============================
check_one_sim <- function(sim_pick = 0) {
  d <- dat %>% filter(sim_id == sim_pick)
  
  sf1 <- survfit(Surv(time, status) ~ 1, data = d)
  na1 <- tibble(time = sf1$time, H_na = -log(sf1$surv))
  
  sp <- survSplit(Surv(time, status) ~ 1, data = d,
                  cut = breaks[-1], episode = "interval",
                  start = "tstart", end = "tstop") %>%
    mutate(interval = factor(interval, labels = cutlabs),
           y = status,
           offset_log_pt = log(pmax(tstop - tstart, .Machine$double.eps)))
  
  fit1 <- glm(y ~ interval - 1 + offset(offset_log_pt),
              data = sp, family = poisson())
  hz <- tibble(interval = levels(sp$interval), h_hat = exp(coef(fit1)))
  
  seg <- tibble(t_left = head(breaks, -1), t_right = tail(breaks, -1)) %>%
    mutate(interval = factor(paste0("(", t_left, ",", t_right, "]"), levels = cutlabs)) %>%
    left_join(hz, by = "interval") %>%
    mutate(width = t_right - t_left,
           H_incr = h_hat * width,
           H_cum_start = c(0, head(cumsum(H_incr), -1)),
           H_cum_end   = cumsum(H_incr))
  
  pw <- bind_rows(
    seg %>% transmute(time = t_left,  H_pw = H_cum_start),
    seg %>% transmute(time = t_right, H_pw = H_cum_end)
  ) %>% arrange(time) %>% distinct(time, .keep_all = TRUE)
  
  p <- ggplot() +
    geom_step(data = na1, aes(time, H_na), linewidth = 0.8) +
    geom_line(data = pw, aes(time, H_pw), linewidth = 0.9, linetype = 2) +
    geom_vline(xintercept = breaks, alpha = 0.4) +
    labs(title = paste0("Cumulative Hazard: Sim ", sim_pick),
         x = "Time (months)", y = "Cumulative hazard") +
    theme_minimal(base_size = 14)
  print(p)
}