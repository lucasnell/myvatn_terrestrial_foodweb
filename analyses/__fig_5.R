

# load packages
suppressPackageStartupMessages({
    library(mtf)
    library(tidyverse)
    library(grid)
})


# This sets plotting device on LAN computer:
if (Sys.getenv("RSTUDIO") == "1" && file.exists(".Rprofile")) {
    source(".Rprofile")
}





dir <- sprintf("~/Box Sync/Iceland Food Web Model/Results/Figures_%s/", Sys.Date())

if (!dir.exists(dir)) dir.create(dir)

#'
#' Note that in `geom_text` below, I added `/ 2.835` to the size arguments to
#' convert from mm to pt.
#'



fig5_caption <- paste("Cumulative midge effect on total bottom-up control",
                      "($BU_{total}$), net top-down control ($TD_{net}$),",
                      "top-down intensification, and top-down alleviation",
                      "as a function of total midge inputs.",
                      "The panels show different combinations of high and",
                      "low predator exploitation and handling times on midges.",
                      "Data are from detritivores; herbivores gave broadly",
                      "similar responses to the cumulative midge effects (Fig. S8).",
                      "Parameter values are as in Table 1, with",
                      "low exploitation $q = 0.8$,",
                      "high exploitation $q = 8$,",
                      "low midge handling time $h_M = 1.3$,",
                      "high midge handling time $h_M = 5.3$,",
                      "pulse duration $w = 20$, and",
                      "pulse rate $b$ ranging from 0.1 to 23 to produce",
                      "different total midge inputs.")



# ================================================================================*
# ================================================================================*

# Simulate data ----

# ================================================================================*
# ================================================================================*



.one_combo <- function(.q, .hM, .b) {

    parlist <- par_estimates %>%
        filter(V==1, H==1, X==1, iI == 10) %>%
        as.list()

    V_gain <- function(V, D, aDV = parlist[["aDV"]]) {
        hD <- parlist[["hD"]]
        (aDV*D*V/(1 + aDV*hD*D)) / V
    }

    V_loss <- function(V, X, H, M, q, hM = parlist[["hM"]], .no_M_to_X = FALSE) {
        aX <- parlist[["aX"]]
        hX <- parlist[["hX"]]
        if (!.no_M_to_X) {
            Vl <- ((aX*V*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M)) / V
        } else {
            Vl <- ((aX*V*X)/(1 + aX*hX*(V + H))) / V
        }
        return(Vl)
    }

    H_gain <- function(P, H, aPH = parlist[["aPH"]]) {
        hP <- parlist[["hP"]]
        (aPH*P*H/(1 + aPH*hP*P)) / H
    }

    H_loss <- function(H, X, V, M, q, hM = parlist[["hM"]], .no_M_to_X = FALSE) {
        aX <- parlist[["aX"]]
        hX <- parlist[["hX"]]
        if (!.no_M_to_X) {
            Hl <- ((aX*H*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M)) / H
        } else {
            Hl <- (aX*H*X)/(1 + aX*hX*(V + H)) / H
        }
        return(Hl)
    }

    # .q = 0.8; .hM = par_estimates$hM[1] * 0.5; .b = 0.1
    # rm(.q, .hM, .b)

    fw <- food_web(tmax = 500, s = 10, b = .b, w = 20,
                   other_pars = list(hM = .hM, q = .q))

    fw <- fw %>%
        spread(pool, N) %>%
        mutate(Vg = V_gain(detritivore, detritus),
               Vl = V_loss(detritivore, predator, herbivore, midge, q = .q, hM = .hM),
               Hg = H_gain(plant, herbivore),
               Hl = H_loss(herbivore, predator, detritivore, midge, q = .q, hM = .hM)) %>%
        summarize(cum_pos_loss_V = sum(Vl[Vl > Vl[1]] - Vl[1]),
                  cum_pos_loss_H = sum(Hl[Hl > Hl[1]] - Hl[1]),
                  cum_neg_loss_V = sum(Vl[Vl < Vl[1]] - Vl[1]),
                  cum_neg_loss_H = sum(Hl[Hl < Hl[1]] - Hl[1]),
                  cum_gain_V = sum(Vg[Vg > Vg[1]] - Vg[1]),
                  cum_gain_H = sum(Hg[Hg > Hg[1]] - Hg[1])) %>%
        mutate(area = 20 * .b,
               q = .q,
               hM = .hM) %>%
        select(area, q, hM, everything()) %>%
        mutate(cum_loss_V = cum_pos_loss_V + cum_neg_loss_V,
               cum_loss_H = cum_pos_loss_H + cum_neg_loss_H,
               # Changing to magnitudes
               cum_neg_loss_V = abs(cum_neg_loss_V),
               cum_neg_loss_H = abs(cum_neg_loss_H)) %>%
        gather("variable", "value", cum_pos_loss_V:cum_loss_H) %>%
        mutate(pool = factor(ifelse(grepl("V$", variable),
                                    "detritivore", "herbivore")),
               direction = factor(ifelse(grepl("loss", variable),
                                         "top-down", "bottom-up")),
               type = case_when(grepl("^cum_pos_loss", variable) ~
                                    "TD intensification",
                                grepl("^cum_neg_loss", variable) ~
                                    "TD alleviation",
                                grepl("^cum_loss_", variable) ~ "TD net",
                                grepl("^cum_gain_", variable) ~ "BU total",
                                TRUE ~ NA_character_) %>%
                   factor(levels = c("BU total", "TD net",
                                     "TD intensification", "TD alleviation"))) %>%
        select(-variable)

    if (sum(is.na(fw$type)) > 0) stop("Error in processing `type` column.")

    return(fw)

}


pulse_df_fig5 <- crossing(.q = c(0.8, 8),
                          .hM = par_estimates$hM[1] * c(0.5, 2),
                          .b = seq(0.1, 23, length.out = 21)) %>%
    pmap_dfr(.one_combo) %>%
    mutate(q = factor(q, levels = sort(unique(q)),
                      labels = sprintf("%s midge\nexploitation",
                                       c("low", "high"))),
           hM = factor(hM, levels = sort(unique(hM)),
                       labels = paste(c("low", "high"),
                                      "midge\nhandling time"))) %>%
    mutate(id = interaction(direction, type))




# ================================================================================*
# ================================================================================*

# Make and write figure ----

# ================================================================================*
# ================================================================================*


fig5_labs <- list(bquote(italic('BU'['total'])),
                  bquote(italic('TD'['net'])),
                  bquote(italic('TD'['intensification'])),
                  bquote(italic('TD'['alleviation'])))


fig5_colors <- c(color_pal()[1], rep("gray60", 3))


td_bu_avail_plot <- pulse_df_fig5 %>%
    filter(pool == "detritivore") %>%
    ggplot(aes(area, value)) +
    geom_hline(yintercept = 0, color = "gray80") +
    geom_line(aes(color = type, linetype = type, group = id), size = 1) +
    geom_text(data = tibble(area =  -75,
                            value = max(pulse_df_fig5$value) * 1.3,
                            hM = pulse_df_fig5$hM %>% unique() %>%
                                sort() %>% rep(each = 2),
                            q = pulse_df_fig5$q %>% unique() %>%
                                sort() %>% rep(2),
                            labs = sprintf("(%s)", letters[1:4])),
              aes(label = labs), hjust = 0, vjust = 1,
              size = 12 / 2.835) +
    facet_grid(hM ~ q) +
    scale_color_manual(NULL, values = fig5_colors, guide = FALSE) +
    scale_linetype_manual(NULL, values = c(1, 1, 3:2),
                          breaks = c("BU total", "TD net",
                                     "TD intensification", "TD alleviation"),
                          labels = fig5_labs) +
    guides(linetype = guide_legend(keywidth = 2, nrow = 2, keyheight = 0.6,
                                   override.aes = list(color = fig5_colors,
                                                       size = 1),
                                   byrow = TRUE)) +
    xlab(expression("Total midge input (" * g ~ N ~ m^{-2} * ")")) +
    ylab("Total midge effect") +
    coord_cartesian(xlim = c(0, 460),
                    ylim = c(-0.404, 4.305),
                    clip = "off") +
    theme(strip.background = element_blank(),
          legend.background = element_blank(),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.justification = "center",
          legend.box = "vertical",
          legend.title = element_text(size = 11, margin = margin(0,0,0,b=-4)),
          legend.text = element_text(size = 10),
          panel.spacing.x = unit(1.5, "lines"),
          panel.spacing.y = unit(2.5, "lines"),
          strip.text.y = element_text(vjust = 1)) +
    NULL


cairo_pdf(filename = paste0(dir, "5-td_bu_avail.pdf"), width = 5, height = 5)
td_bu_avail_plot
dev.off()






# ================================================================================*
# ================================================================================*

# Write caption ----

# ================================================================================*
# ================================================================================*



cf <- file(paste0(dir, "captions.txt"), "a")

writeLines(paste0("Figure 5. ", fig5_caption, "\n"), cf)

close(cf)
