
# Run the simulations analyzed in long-term_dynamics.R


# load packages
suppressPackageStartupMessages({
    library(mtf)
    library(tidyverse)
    library(parallel)
})


# Number of cores to use:
n_cores <- as.integer(detectCores() * 0.9)

# ------------------------
# Summary functions
# ------------------------

# Time to max:
to_max <- function(x, s, w) {
    x <- x[s:length(x)]
    which(x == max(x))[1] - w
}
# Time to min:
to_min <- function(x, s, w) {
    x <- x[s:length(x)]
    which(x == min(x))[1] - w
}
# Cumulative N:
cum_N <- function(x) sum(x)
# Min N below equilibrium
min_below <- function(x) min(x - x[1])
# Max N above equilibrium
max_above <- function(x) max(x - x[1])
# Return time
return_time <- function(x, s, w, thresh = 1e-1, ...) {
    n <- length(x)
    x_dr <- diff(range(x))
    y <- abs(rev((x[(s+w):n] - x[1]) / x_dr))
    ind <- which(y > thresh)
    if (length(ind) == 0) return(NA_real_)
    if (ind[1] == length(y)) return(NaN)
    ind <- ind[1] - 1
    ind <- n - ind + 1
    return(ind - (s + w))
}


# ------------------------
# Simulation function
# ------------------------

# Runs simulations for one combination of parameters
one_combo <- function(row_i) {
    .w <- row_i$w
    .b <- row_i$b
    .lM <- row_i$lM
    # .mM <- row_i$mM
    .mM <- 0.5
    .aM <- row_i$aM
    fw <- food_web(tmax = 250, s = 10, b = .b, w = .w, .lM = .lM, .aM = .aM, .mM = .mM)
    fw <- fw %>%
        filter(pool != "midge") %>%
        mutate(pool = droplevels(pool)) %>%
        group_by(pool) %>%
        summarize(to_max = to_max(N, s = 10, w = .w),
                  to_min = to_min(N, s = 10, w = .w),
                  cum_N = cum_N(N),
                  min_below = min_below(N),
                  max_above = max_above(N),
                  return_time = return_time(N, s = 10, w = .w)) %>%
        ungroup() %>%
        mutate(w = .w, b = .b, lM = .lM, mM = .mM, aM = .aM, area = b * w) %>%
        select(w, b, lM, area, everything())
    return(fw)
}


# ------------------------
# Combinations of parameter values
# ------------------------

pulse_pars <- expand.grid(w = seq(10, 30, length.out = 100),
                          b = seq(0.1, 100, length.out = 100),
                          lM = c(0.1, 0.325, 0.55),
                          # mM = c(0.1, 0.5, 2.5),
                          aM = c(1e-3, 5e-3, 1e-2)) %>%
    split(row(.)[,1])


# ------------------------
# Run and summarize simulations
# ------------------------

pulse_df <- mclapply(pulse_pars, one_combo, mc.cores = n_cores)
pulse_df <- bind_rows(pulse_df) %>%
    select(w, b, lM, mM, aM, pool, everything()) %>%
    mutate(pool = fct_recode(pool, soil = "nitrogen"))


# ------------------------
# Make sure the following return TRUE:
# ------------------------

# Verify that there aren't any NaN values (i.e., that they all returned to equilibrium):
pulse_df %>%
    gather("param", "value", to_max:return_time, factor_key = TRUE) %>%
    filter(is.nan(value)) %>%
    nrow() %>%
    `==`(0)
# Verify that there aren't any NA values (i.e., that they all were perturbed enough):
pulse_df %>%
    gather("param", "value", to_max:return_time, factor_key = TRUE) %>%
    filter(is.na(value)) %>%
    nrow() %>%
    `==`(0)
# If -to_min or -to_max is > w, that means that those functions malfunctioned
pulse_df %>%
    filter(-to_min > w || -to_max > w) %>%
    nrow() %>%
    `==`(0)


# ------------------------
# Make sure the following values generally make sense:
# ------------------------

# Summary of different parameters:
pulse_df %>%
    gather("param", "value", to_max:return_time, factor_key = TRUE) %>%
    group_by(param) %>%
    summarize(min = min(value, na.rm = TRUE), max = max(value, na.rm = TRUE))


write_csv(pulse_df, "data-raw/pulse_data.csv")

