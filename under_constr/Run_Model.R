# ======================================
# Load Package
# ======================================

library(mtf)
library(tidyverse)



# # Test midge pulse:
# plot(test_midges(1000, a=1e9, b=10, r=400, w=40, d=1), type = 'l', ylab = "iM")


web_output <- food_web(tmax = 1000, a = 1e9, b = 0, r = 400, w = 40, d = 1)


# Plot absolute N
web_output %>%
    group_by(pool) %>%
    # Define 'minb' to set the minimum value for the y-axis.
    # This allows different y-scales for different facets, with the ymin set to 0
    mutate(minb = 0) %>%
    ggplot(aes(time, N)) +
    facet_wrap(~pool, scales="free_y") +
    # The horizontal lines show the initial states
    geom_hline(data = web_output %>% filter(time == min(time)),
               aes(yintercept=N), color="firebrick") +
    geom_line(size = 1) +
    geom_point(aes(time, minb), shape="") +
    theme_classic()



modify_b <- map_dfr(c(0, 10, 20, 40),
                    function(.b) {
                        food_web(tmax = 1000, a = 1e9, b = .b, r = 400, w = 40, d = 1) %>%
                            mutate(b = .b)})

modify_b %>%
    mutate(b = factor(b)) %>%
    filter(pool != "midge") %>%
    group_by(pool, b) %>%
    # Scale the N relative to the initial state
    mutate(N_scale = N/N[1]) %>%
    ggplot(aes(time, N_scale, color = b)) +
    facet_wrap(~pool, scales = 'free_y') +
    # geom_hline(yintercept = 1, color="firebrick4") +
    geom_line(size = 1) +
    theme_classic() +
    scale_x_continuous(breaks = c(0, 500, 1000))



# Relative to Equilibrium
# Note that this gives weird results when there is no deviation from equilibrium
# This is probably due to small numerical errors, but is not a big issue
web_output %>%
    filter(pool!="midge") %>%
    group_by(pool) %>%
    # Scale the N relative to the initial state
    mutate(N_scale = N/N[1]) %>%
    ungroup() %>%
    ggplot(aes(time, N_scale)) +
    facet_wrap(~pool) +
    geom_hline(yintercept = 1, color="firebrick4") +
    geom_line(size = 1) +
    theme_classic()






