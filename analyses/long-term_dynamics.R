# ======================================
# Preliminaries
# ======================================

# load packages
suppressPackageStartupMessages({
    library(mtf)
    library(tidyverse)
    library(forcats)
})

pulse_df <- read_csv("data-raw/pulse_data.csv",
                     col_types = cols(
                         w = "f",
                         b = "f",
                         lM = "f",
                         pool = "f",
                         mM = "f",
                         aM = "f",
                         .default = "d"))
pulse_df <- pulse_df %>%
    select(w, b, lM, mM, aM, pool, everything()) %>%
    mutate(pool = fct_recode(pool, soil = "nitrogen"))


# ------------------------
# Plotting functions
# ------------------------

# # Construct a heatmap plot
# heat_plot <- function(.df, .col, .pool) {
#     .col <- enquo(.col)
#     .df <- .df %>%
#         filter(pool == .pool) %>%
#         mutate_at(vars(lM, mM), factor)
#     midpt <- .df %>% select(!!.col) %>% .[[1]] %>% range(na.rm = TRUE) %>% median()
#     color_lab_ <- gsub("\\_", " ", gsub("~", "", deparse(.col))) %>%
#         trimws() %>% tools::toTitleCase()
#     .df %>%
#         ggplot(aes(w, b)) +
#         geom_tile(aes(fill = !!.col)) +
#         geom_line(data = tibble(w = rep(seq(10, 30, length.out = 100), 5),
#                                 a = rep(seq(250, 1250, 250), each = 100),
#                                 b = a / w),
#                   aes(group = factor(a)), size = 1) +
#         scale_fill_gradient2(color_lab_,
#                              low = "dodgerblue", mid = "white", high = "firebrick",
#                              midpoint = midpt) +
#         facet_grid(mM ~ lM, labeller = label_both) +
#         ggtitle(.pool) +
#         coord_cartesian(ylim = c(0, 100), xlim = c(10, 30)) +
#         ylab("Rate of input") +
#         xlab("Duration of input")
# }
# Plot of area vs another variable
area_plot <- function(.df, .col, .pools, .ylab = NULL) {
    .col <- enquo(.col)
    .df <- .df %>%
        filter(pool %in% .pools, mM == 0.5, lM == 0.5) %>%
        # filter(pool %in% .pools, mM == 0.5, aM == 1e-2) %>%
        mutate(w = as.numeric(paste(w)), pool = factor(paste(pool), levels = .pools))
    if (is.null(.ylab)) {
        .ylab <- gsub("\\_", " ", gsub("~", "", deparse(.col))) %>%
            trimws() %>% tools::toTitleCase()
    }
    .df %>%
        ggplot(aes(area, return_time)) +
        geom_hline(yintercept = 75, linetype = 2, color = "gray70") +
        geom_point(aes(color = w), shape = 1, alpha = 0.5, na.rm = TRUE) +
        scale_color_gradient2("Width",
                              low = "firebrick", mid = "gray80", high = "dodgerblue",
                              midpoint = 20) +
        # ggtitle(.pool) +
        facet_grid(pool ~ aM) +
        # facet_grid(pool ~ lM) +
        ylab(.ylab) +
        xlab("Area under curve")
}



# ------------------------
# Plots
# ------------------------



area_plot(pulse_df, return_time, c("detritivore", "herbivore", "predator"))
area_plot(pulse_df, return_time, c("plant", "soil", "detritus"))





heat_plot(pulse_df, to_max, "detritivore")
heat_plot(pulse_df, to_min, "detritivore")
heat_plot(pulse_df, cum_N, "detritivore")



pulse_df %>%
    filter(pool == "detritivore", lM == 0.1) %>%
    # mutate(lM = factor(lM)) %>%
    # ggplot(aes(w, return_time)) +
    ggplot(aes(b, return_time)) +
    # geom_point(aes(color = w)) +
    # geom_line(aes(color = b, group = factor(b)), size = 1) +
    geom_line(aes(color = w, group = factor(w)), size = 1) +
    # facet_grid(~ lM) +
    # facet_wrap(~ w, nrow = 5) +
    scale_color_gradient2("Width",
                          low = "firebrick", mid = "gray80", high = "dodgerblue",
                          midpoint = 20) +
    # scale_color_gradient2("Intensity",
    #                       low = "firebrick", mid = "gray80", high = "dodgerblue",
    #                       midpoint = 50) +
    NULL


pulse_df %>%
    filter(pool == "detritivore") %>%
    mutate_at(vars(lM, mM), factor) %>%
    ggplot(aes(cum_N, max_above)) +
    geom_line(aes(color = w, group = factor(b)), size = 1) +
    facet_grid(mM ~ lM, labeller = label_both) +
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
    # geom_vline(data = tibble(xint = c(200, 900, 1100), lM = factor(c(0.1, 0.5, 0.9))),
    #           aes(xintercept = xint), linetype = 2)
area_plot(pulse_df, return_time, "predator")

area_plot(pulse_df, return_time, "nitrogen")
area_plot(pulse_df, return_time, "detritus")
area_plot(pulse_df, return_time, "plant")

