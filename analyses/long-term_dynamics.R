# ======================================
# Preliminaries
# ======================================

# load packages
library(mtf)
library(tidyverse)
library(parallel)


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
    y <- abs(rev((x[(s+w):n] - x[1]) / x[1]))
    ind <- which(y > thresh)
    if (length(ind) == 0) return(NA_real_)
    if (ind[1] == length(y)) return(NaN)
    ind <- ind[1] - 1
    ind <- n - ind + 1
    return(ind - (s + w))
}
# Runs simulations for one combination of parameters
one_combo <- function(row_i) {
    .w <- row_i$w
    .b <- row_i$b
    .aM <- row_i$aM
    fw <- food_web(tmax = 250, s = 10, b = .b, w = .w, .lM = 0.1, .aM = .aM)
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
# Construct a heatmap plot
heat_plot <- function(.df, .col, .pool) {
    .col <- enquo(.col)
    .df <- .df %>%
        filter(pool == .pool) %>%
        mutate(aM = factor(aM))
    midpt <- .df %>% select(!!.col) %>% .[[1]] %>% range(na.rm = TRUE) %>% median()
    color_lab_ <- gsub("\\_", " ", gsub("~", "", deparse(.col))) %>%
        trimws() %>% tools::toTitleCase()
    .df %>%
        ggplot(aes(w, b)) +
        geom_tile(aes(fill = !!.col)) +
        geom_line(data = tibble(w = rep(seq(10, 30, length.out = 100), 5),
                                a = rep(seq(250, 1250, 250), each = 100),
                                b = a / w),
                  aes(group = factor(a)), size = 1) +
        scale_fill_gradient2(color_lab_,
                             low = "dodgerblue", mid = "white", high = "firebrick",
                             midpoint = midpt) +
        facet_grid(~ aM) +
        ggtitle(.pool) +
        coord_cartesian(ylim = c(0, 100), xlim = c(10, 30)) +
        ylab("Rate of input") +
        xlab("Duration of input")
}
# Plot of area vs another variable
area_plot <- function(.df, .col, .pool) {
    .col <- enquo(.col)
    .df <- .df %>%
        filter(pool == .pool) %>%
        mutate(aM = factor(aM))
    ylab_ <- gsub("\\_", " ", gsub("~", "", deparse(.col))) %>%
        trimws() %>% tools::toTitleCase()
    .df %>%
        ggplot(aes(area, !!.col)) +
        geom_point(aes(color = w), shape = 1, alpha = 0.5) +
        scale_color_gradient2("Width",
                              low = "firebrick", mid = "gray80", high = "dodgerblue",
                              midpoint = 20) +
        facet_grid(~ aM) +
        ggtitle(.pool) +
        ylab(ylab_) +
        xlab("Area under curve")
}

# ------------------------
# Combinations of parameter values
# ------------------------

pulse_pars <- expand.grid(w = seq(10, 30, length.out = 10),
                          b = seq(0.1, 100, length.out = 10),
                          aM = c(1e-4, 1e-2, 1)) %>%
    split(row(.)[,1])


# ------------------------
# Run and summarize simulations
# ------------------------

# Takes ~45 min w/ 3 cores
pulse_df <- mclapply(pulse_pars, one_combo, mc.cores = 3)
pulse_df <- bind_rows(pulse_df)

# write_csv(pulse_df, "data-raw/pulse_data.csv")

# pulse_df <- read_csv("data-raw/pulse_data.csv", col_types = "ddddfdddddd")

# Verify that there aren't any NaN values (i.e., that they all returned to equilibrium):
pulse_df %>%
    gather("param", "value", to_max:return_time, factor_key = TRUE) %>%
    filter(is.nan(value)) %>%
    nrow() %>%
    `==`(0)
# If -to_min or -to_max is > w, that means that those functions malfunctioned
pulse_df %>%
    filter(-to_min > w || -to_max > w) %>%
    nrow() %>%
    `==`(0)
# Summary of different parameters:
pulse_df %>%
    gather("param", "value", to_max:return_time, factor_key = TRUE) %>%
    group_by(param) %>%
    summarize(min = min(value, na.rm = TRUE), max = max(value, na.rm = TRUE))



# ------------------------
# Plots
# ------------------------



heat_plot(pulse_df, to_max, "detritivore")
heat_plot(pulse_df, to_min, "detritivore")
heat_plot(pulse_df, cum_N, "detritivore")


pulse_df %>%
    filter(pool == "detritivore", aM == 1, w == 22) %>%
    select(b, return_time)

pulse_df %>%
    filter(pool == "detritivore", aM == 1) %>%
    # mutate(aM = factor(aM)) %>%
    # ggplot(aes(w, return_time)) +
    ggplot(aes(b, return_time)) +
    # geom_point(aes(color = w)) +
    # geom_line(aes(color = b, group = factor(b)), size = 1) +
    geom_line(aes(color = w, group = factor(w)), size = 1) +
    # facet_grid(~ aM) +
    facet_wrap(~ w, nrow = 5) +
    scale_color_gradient2("Width",
                          low = "firebrick", mid = "gray80", high = "dodgerblue",
                          midpoint = 20) +
    # scale_color_gradient2("Intensity",
    #                       low = "firebrick", mid = "gray80", high = "dodgerblue",
    #                       midpoint = 50) +
    NULL


pulse_df %>%
    filter(pool == "detritivore") %>%
    mutate(aM = factor(aM)) %>%
    ggplot(aes(cum_N, max_above)) +
    geom_line(aes(color = w, group = factor(b)), size = 1) +
    facet_grid(~ aM) +
    scale_color_gradient2("Width",
                          low = "firebrick", mid = "gray80", high = "dodgerblue",
                          midpoint = 20)


heat_plot(pulse_df, min_below, "detritivore")
heat_plot(pulse_df, max_above, "detritivore")

heat_plot(pulse_df, return_time, "detritivore")


heat_plot(pulse_df, return_time, "herbivore")
heat_plot(pulse_df, return_time, "predator")






area_plot(pulse_df, return_time, "detritivore") +
    geom_vline(xintercept = 350, linetype = 2)
area_plot(pulse_df, return_time, "herbivore")  # +
    # geom_vline(data = tibble(xint = c(200, 900, 1100), aM = factor(c(1e-4, 0.01, 1))),
    #           aes(xintercept = xint), linetype = 2)
area_plot(pulse_df, return_time, "predator")

area_plot(pulse_df, return_time, "nitrogen")
area_plot(pulse_df, return_time, "detritus")
area_plot(pulse_df, return_time, "plant")

