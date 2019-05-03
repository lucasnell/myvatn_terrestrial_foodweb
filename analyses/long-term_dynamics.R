

# load packages
suppressPackageStartupMessages({
    library(mtf)
    library(tidyverse)
    library(grid)
    library(pbmcapply)
})




#'
#' Note that in `geom_text` below, I added `/ 2.835` to the size arguments to
#' convert from mm to pt.
#'




full_pulse_df <- read_csv("data-raw/pulse_data.csv",
                         col_types = cols(
                             w = "f",
                             pool = "f",
                             .default = "d")) %>%
        mutate_at(vars(b, f), ~ factor(signif(., digits = 6)))

make_vec <- function(cname) rep(signif(par_combs[[cname]], digits = 6), each = 6)

par_combs <- expand.grid(w = seq(10, 25, length.out = 25),
                         b = seq(0.1, 40, length.out = 25),
                         f = seq(8e-3, 8e-2, length.out = 25),
                         mM = par_estimates$mM[1] * c(0.5, 1, 2),
                         hM = par_estimates$hM[1] * c(0.5, 1, 2),
                         # Plant / herbivore uptake rates:
                         aDV = c(par_estimates$aDV[1], par_estimates$aPH[1]),
                         aPH = c(par_estimates$aPH[1], par_estimates$aDV[1])) %>%
    mutate_at(vars(b, f, mM, hM, aDV, aPH), ~ signif(., digits = 6))

pulse_df <- full_pulse_df %>%
    mutate(mM = make_vec("mM"),
           hM = make_vec("hM"),
           aDV = make_vec("aDV"),
           aPH = make_vec("aPH")) %>%
    filter(mM == signif(par_estimates$mM[1], digits = 6),
           hM == signif(par_estimates$hM[1], digits = 6),
           aDV == signif(par_estimates$aDV[1], digits = 6),
           aPH == signif(par_estimates$aPH[1], digits = 6)) %>%
    select(-mM, -hM, -aDV, -aPH)



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
        mutate(f = f %>% paste() %>% as.numeric())
    if (flip_sign) {
        .df <- .df %>% mutate_at(vars(!!.col), ~ . * -1)
    }
    p <- .df %>%
        ggplot(aes(f, area)) +
        geom_tile(aes(fill = !!.col)) +
        geom_contour(aes(z = !!.col), color = "white", bins = 7) +
        scale_fill_viridis_c(.legend_title, option = .pal_opt, breaks = .legend_breaks,
                             labels = .legend_labels) +
        # ggtitle(.title) +
        ylab("Total midge input") +
        scale_x_continuous("Midge availability", breaks = c(0.04, 0.08)) +
        theme(legend.title = element_text(size = 10),
              legend.key.width = grid::unit(0.02, "npc"),
              plot.title = element_text(hjust = 0, size = 12))
    if (!is.null(.facet_lab)) {
        p <- p  +
            geom_text(data = tibble(f = 0.008, area = 750, lab = .facet_lab),
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





hm1 <- heat_plot(cum_pos_loss_V, .pal_opt = "B", .facet_lab = "a",
                 .legend_breaks = c(0.1, 1.5),
                 keep_y = FALSE, keep_x = FALSE,
                 .legend_title = "Top-down\nintensification")
hm2 <- heat_plot(cum_neg_loss_V, .pal_opt = "D", flip_sign = TRUE,  .facet_lab = "b",
                 .legend_breaks = c(0.01, 0.075),
                 .legend_title = "Top-down\nalleviation", keep_x = FALSE)
hm3 <- heat_plot(cum_gain_V, .legend_title = "Bottom-up",  .facet_lab = "c",
                 keep_y = FALSE,
                 .pal_opt = "C") +
    scale_fill_gradient("Bottom-up", high = "deepskyblue", low = "black",
                        breaks = c(1, 5))
# heat_plot(cum_gain_H, "herbivore", .title = "Enhancement of bottom-up effect", flip_color = TRUE)


hmg1 <- ggplotGrob(hm1)
hmg2 <- ggplotGrob(hm2)
hmg3 <- ggplotGrob(hm3)

hm <- rbind(hmg1, hmg2, hmg3, size = "first")
hm$widths <- unit.pmax(hmg1$widths, hmg2$widths, hmg3$widths)



# grid.newpage()
# grid.draw(hm)


# pdf(file = "~/Desktop/5-heatmaps.pdf", width = 4.5, height = 9)
# grid.newpage()
# grid.draw(hm)
# dev.off()




td_bu_plot <- pulse_df %>%
    mutate_at(vars(f, b), function(x) x %>% paste() %>% as.numeric()) %>%
    filter(w == 20,
           pool %in% c("detritivore", "herbivore"),
           f %in% as.numeric(quantile(f, 0.5))) %>%
    group_by(f, b, pool) %>%
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
    geom_text(data = tibble(cum_gain = c(6, 2.5),
                            cum_loss = rep(1.4, 2),
                            pool = factor(c("detritivore", "herbivore"),
                                          levels = c("detritivore", "herbivore"))),
              aes(label = pool, color = pool), hjust = 0, vjust = 1, size = 10 / 2.835) +
    scale_color_manual(values = color_pal()[1:2], guide = FALSE) +
    scale_size_continuous("Midge pulse\nintensity", range = c(0.5, 6),
                          breaks = c(1, 5, 25)) +
    xlab("Bottom-up effect") +
    ylab("Total top-down effect") +
    theme(legend.position = c(0.85, 0.25), legend.background = element_blank()) +
    NULL





# ggsave(filename = "~/Desktop/6-top_vs_bottom.pdf", td_bu_plot, width = 7, height = 4)


