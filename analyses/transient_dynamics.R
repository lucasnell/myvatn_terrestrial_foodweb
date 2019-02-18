# ======================================
# Preliminaries
# ======================================

# load packages
library(mtf)
library(tidyverse)




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
    y <- abs(rev((x[(s+w):n] - x[1]) / x[1]))
    ind <- tail(which(y <= thresh), 1)
    if (length(ind) == 0) return(NA_real_)
    ind <- n - ind + 1
    return(ind - (s + w))
}




one_combo <- function(row_i) {
    .w <- row_i$w
    .b <- row_i$b
    .aM <- row_i$aM
    fw <- food_web(tmax = 200, s = 10, b = .b, w = .w, .lM = 0.1, .aM = .aM)
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
        mutate(w = .w, b = .b, aM = .aM, area = b * w) %>%
        select(w, b, aM, area, everything())
    return(fw)
}

pulse_pars <- expand.grid(w = seq(10, 30, length.out = 100),
                          b = seq(0.1, 100, length.out = 100),
                          aM = c(1e-4, 1e-2, 1)) %>%
    split(row(.)[,1])

library(parallel)
t0 <- Sys.time()
pulse_df <- mclapply(pulse_pars, one_combo, mc.cores = 3)
Sys.time() - t0; rm(t0)



heat_plot <- function(.df, .col, .pool) {
    .col <- enquo(.col)
    .df <- .df %>%
        filter(pool == .pool) %>%
        mutate(aM = factor(aM))
    midpt <- .df %>% select(!!.col) %>% .[[1]] %>% range() %>% median()
    .df %>%
        ggplot(aes(w, b)) +
        geom_tile(aes(fill = !!.col)) +
        geom_line(data = tibble(w = rep(seq(10, 30, length.out = 100), 5),
                                a = rep(seq(250, 1250, 250), each = 100),
                                b = a / w),
                  aes(group = factor(a)), size = 1) +
        scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick",
                             midpoint = midpt) +
        facet_grid(~ aM) +
        coord_cartesian(ylim = c(0, 100), xlim = c(10, 30))
}

heat_plot(pulse_pars, return_time, "detritivore")

