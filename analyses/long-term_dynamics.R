
# load packages
suppressPackageStartupMessages({
    library(mtf)
    library(tidyverse)
    library(forcats)
    library(grid)
    library(parallel)
})


#'
#' Note that in `geom_text` below, I added `/ 2.835` to the size arguments to
#' convert from mm to pt.
#'





pulse_df <- read_csv("data-raw/pulse_data.csv",
                     col_types = cols(
                         w = "f",
                         b = "f",
                         aM = "f",
                         pool = "f",
                         .default = "d"))

# < order of colors: green, red, purple, pink, light green, yellow >
# rows correspond to `RColorBrewer::brewer.pal(6, "Dark2")`
rgb_mat <- rbind(c(27,158,119), c(217,95,2), c(117,112,179),
                 c(231,41,138), c(102,166,30), c(230,171,2))
# Switching detritus with plant and detritivore with herbivore
# Order is now "detritivore", "herbivore", "predator", "soil", "detritus", "plant"
rgb_mat <- rgb_mat[c(2, 5, 3, 6, 4, 1),]
# Multipliers for each color, < 1 makes it darker
rgb_mults <- c(0.5, 1.3, 1.1,
               1.05, 0.9, 1)
color_pal <- apply(rgb_mat * matrix(rgb_mults, 6, 3), 1,
                   function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))




heat_plot <- function(.col, .facet_lab = NULL,
                      .pool = "detritivore", .df = NULL, .w = NULL, .title = NULL,
                      .legend_title = NULL,
                      .legend_breaks = waiver(), .legend_labels = waiver(),
                      flip_sign = FALSE, flip_color = FALSE,
                      keep_x = TRUE, keep_y = TRUE, .pal_opt = "C") {
    if (is.null(.df)) .df <- pulse_df
    if (is.null(.w)) .w <- .df$w %>% paste() %>% as.numeric() %>% median()
    .col <- enquo(.col)
    # if (is.null(.title)) .title <- .pool
    .df <- .df %>%
        filter(pool == .pool, w == .w) %>%
        mutate(aM = aM %>% paste() %>% as.numeric())
    if (flip_sign) {
        .df <- .df %>% mutate_at(vars(!!.col), ~ . * -1)
    }
    p <- .df %>%
        ggplot(aes(aM, area)) +
        geom_tile(aes(fill = !!.col)) +
        geom_contour(aes(z = !!.col), color = "white", bins = 7) +
        scale_fill_viridis_c(.legend_title, option = .pal_opt, breaks = .legend_breaks,
                             labels = .legend_labels) +
        # ggtitle(.title) +
        ylab("Total midge input") +
        scale_x_continuous("Attack rate on midges", breaks = c(0.005, 0.01)) +
        theme(legend.title = element_text(size = 10),
              legend.key.width = grid::unit(0.02, "npc"),
              plot.title = element_text(hjust = 0, size = 12))
    if (!is.null(.facet_lab)) {
        p <- p +
            geom_text(data = tibble(aM = 0.001, area = 2100, lab = .facet_lab),
                      aes(label = lab), hjust = 0, vjust = 0, size = 12 / 2.835)
    }
    if (!keep_x) {
        p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    }
    if (!keep_y) {
        p <- p + theme(axis.title.y = element_blank())
    }
    return(p)
}



#'
#' - Overall label for rows in first two plots?
#' - Labeling herbivores in herbivore-only plot?
#'
#'

hm1 <- heat_plot(cum_pos_loss_V, .pal_opt = "B", .facet_lab = "a",
                 .legend_breaks = c(0.1, 2),
                 keep_y = FALSE,
                 .legend_title = "Top-down\nintensification", keep_x = FALSE)
hm2 <- heat_plot(cum_neg_loss_V, .pal_opt = "D", flip_sign = TRUE,  .facet_lab = "b",
                 .legend_breaks = c(0.01, 0.2),
                 .legend_title = "Top-down\nalleviation", keep_x = FALSE)
hm3 <- heat_plot(cum_gain_V, .legend_title = "Bottom-up",  .facet_lab = "c",
                 keep_y = FALSE,
                 .pal_opt = "C") +
    scale_fill_gradient("Bottom-up", high = "deepskyblue", low = "black",
                        breaks = c(1, 9))
# heat_plot(cum_gain_H, "herbivore", .title = "Enhancement of bottom-up effect", flip_color = TRUE)

hmg1 <- ggplotGrob(hm1)
hmg2 <- ggplotGrob(hm2)
hmg3 <- ggplotGrob(hm3)

hm <- rbind(hmg1, hmg2, hmg3, size = "first")
hm$widths <- unit.pmax(hmg1$widths, hmg2$widths, hmg3$widths)


grid.newpage()
grid.draw(hm)


# pdf(file = "~/Desktop/5-heatmaps.pdf", width = 4.5, height = 9)
# grid.newpage()
# grid.draw(hm)
# dev.off()



parlist <- par_estimates %>%
    filter(V==1, H==1, R==1, iN == 10) %>%
    as.list()
V_loss <- function(V, R, H, M, aM) {
    aR <- parlist[["aR"]]
    hVH <- parlist[["hVH"]]
    hM <- parlist[["hM"]]
    ((aR*V*R)/(1 + aR*hVH*(V + H) + aM*hM*M)) / V
}


# Takes ~45 sec w/ 4 cores
td_input_sims <- crossing(aM = c(0.00175, 0.00550, 0.00925),
                          w = c(10, 30),
                          area = seq(0, 3000, length.out = 100)) %>%
    split(1:nrow(.)) %>%
    mclapply(
        function(.df) {
            .aM <- .df$aM
            .w <- .df$w
            .area <- .df$area
            .b <- .area / .w
            fw <- food_web(tmax = 250, s = 10, b = .b, w = .w,
                           other_pars = list(aM = .aM))
            Vl <- fw %>%
                spread(pool, N) %>%
                mutate(Vl = V_loss(detritivore, predator, herbivore, midge, .aM)) %>%
                .[["Vl"]]
            .td <- sum(Vl[Vl > Vl[1]] - Vl[1]) + sum(Vl[Vl < Vl[1]] - Vl[1])
            return(mutate(.df, td = .td))
        }, mc.cores = detectCores()) %>%
    bind_rows()


td_input_plot <- td_input_sims %>%
    mutate_at(vars(aM, area), function(x) x %>% paste() %>% as.numeric()) %>%
    mutate(aM = factor(aM, levels = sort(unique(aM)),
                       labels = c("low", "mid", "high")),
           b = area / w,
           w = factor(w)) %>%
    ggplot(aes(area)) +
    geom_line(aes(y = td, color = aM, linetype = w), size = 1) + # net top-down
    # Manually specify colors from viridis's 5-color "inferno" palette
    scale_color_manual("Attack rate",
                       values = c("#000004FF", "#4C0C6BFF", "#A82E5FFF",
                                  "#EF6C23FF", "#F6D645FF")[c(1, 4, 5)]) +
    scale_linetype_manual("Pulse width", values = c(2, 1)) +
    xlab("Total midge input") +
    ylab("Total top-down effect") +
    theme(legend.position = c(0.75, 0.25), legend.box = "horizontal",
          legend.title = element_text(size = 10),
          legend.key.width = grid::unit(0.075, "npc"))


# ggsave(filename = "~/Desktop/6-td_input.pdf", td_input_plot, width = 7, height = 4)





td_bu_plot <- pulse_df %>%
    mutate_at(vars(aM, b), function(x) x %>% paste() %>% as.numeric()) %>%
    filter(w == 20,
           pool %in% c("detritivore", "herbivore"),
           aM %in% as.numeric(quantile(aM, 0.5))) %>%
    group_by(aM, b, pool) %>%
    summarize(cum_pos_loss = ifelse(pool[1] == "detritivore", cum_pos_loss_V[1],
                                    cum_pos_loss_H[1]),
              cum_neg_loss = ifelse(pool[1] == "detritivore", cum_neg_loss_V[1],
                                    cum_neg_loss_H[1]),
              cum_loss = cum_pos_loss + cum_neg_loss,
              cum_gain = ifelse(pool[1] == "detritivore", cum_gain_V[1],
                                cum_gain_H[1])) %>%
    ungroup() %>%
    ggplot(aes(cum_gain, cum_loss)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray50") +
    geom_line(aes(color = pool), size = 1) +
    geom_point(aes(size = b, color = pool)) +
    geom_text(data = tibble(cum_gain = c(9.75, 4.25),
                            cum_loss = rep(2.45, 2),
                            pool = factor(c("detritivore", "herbivore"),
                                          levels = c("detritivore", "herbivore"))),
              aes(label = pool, color = pool), hjust = 1, vjust = 1, size = 10 / 2.835) +
    scale_color_manual(values = color_pal[1:2], guide = FALSE) +
    scale_size_continuous("Midge pulse\nintensity", range = c(0.5, 6),
                          breaks = c(10, 100)) +
    xlab("Bottom-up effect") +
    ylab("Total top-down effect") +
    theme(legend.position = c(0.85, 0.25), legend.background = element_blank()) +
    NULL





# ggsave(filename = "~/Desktop/7-top_vs_bottom.pdf", td_bu_plot, width = 7, height = 4)


