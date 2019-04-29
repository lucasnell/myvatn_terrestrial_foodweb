
# Run the simulations analyzed in long-term_dynamics.R


if (!require(pbmcapply)) install.packages("pbmcapply")


# load packages
suppressPackageStartupMessages({
    library(mtf)
    library(tidyverse)
    library(pbmcapply)
})

# Number of cores to use:
n_cores <- as.integer(detectCores() * 0.9)

# ------------------------
# Summary functions
# ------------------------

parlist <- par_estimates %>%
        filter(V==1, H==1, R==1, iN == 10) %>%
        as.list()
V_gain <- function(V, D, aDV) {
    hD <- parlist[["hD"]]
    (aDV*D*V/(1 + aDV*hD*D)) / V
}
V_loss <- function(V, R, H, M, f, hM) {
    aR <- parlist[["aR"]]
    hVH <- parlist[["hVH"]]
    ((aR*V*R)/(1 + aR*hVH*(V + H) + (aR * f)*hM*M)) / V
}
H_gain <- function(P, H, aPH) {
    hP <- parlist[["hP"]]
    (aPH*P*H/(1 + aPH*hP*P)) / H
}
H_loss <- function(H, R, V, M, f, hM) {
    aR <- parlist[["aR"]]
    hVH <- parlist[["hVH"]]
    ((aR*H*R)/(1 + aR*hVH*(V + H) + (aR * f)*hM*M)) / H
}

M_flux_R <- function(H, R, V, M, f, hM) {
    aR <- parlist[["aR"]]
    hVH <- parlist[["hVH"]]
    ((aR * f)*M*R)/(1 + aR*hVH*(V + H) + (aR * f)*hM*M) / M
}
M_flux_D <- function(M, mM) {
    mM*M
}








# ------------------------
# Simulation function
# ------------------------

# Runs simulations for one combination of parameters
one_combo <- function(row_i) {
    .w <- row_i$w
    .b <- row_i$b
    .other_pars <- as.list(unlist(row_i))
    .other_pars$w <- NULL
    .other_pars$b <- NULL

    # Changing aDV or aPH changed equil. pool sizes, so need to find equil_pools objects
    # that specify the new sizes if either arg is not the default:
    if (.other_pars$aDV == par_estimates$aDV[1] &&
        .other_pars$aPH == par_estimates$aPH[1]) {  # defaults
        .ep_obj <- NULL
    } else {
        .ep_obj <- ep_df %>%
            filter(aDV == .other_pars$aDV, aPH == .other_pars$aPH) %>%
            .[["pools"]]
        if (length(.ep_obj) == 0) {
            stop("\nDesired combination of aDV and aPH isn't available")
        }
        .ep_obj <- .ep_obj[[1]]
    }

    fw <- food_web(tmax = 250, s = 10, b = .b, w = .w, ep_obj = .ep_obj,
                   other_pars = .other_pars)

    fw <- fw %>%
        spread(pool, N) %>%
        mutate(Vg = V_gain(detritivore, detritus, .other_pars$aDV),
               Vl = V_loss(detritivore, predator, herbivore, midge, .other_pars$f,
                           .other_pars$hM),
               Hg = H_gain(plant, herbivore, .other_pars$aPH),
               Hl = H_loss(herbivore, predator, detritivore, midge, .other_pars$f,
                           .other_pars$hM),
               MfR = M_flux_R(herbivore, predator, detritivore, midge, .other_pars$f,
                        .other_pars$hM),
               MfD = M_flux_D(midge, .other_pars$mM)) %>%
        gather("pool", "N", soil:midge) %>%
        filter(pool != "midge") %>%
        mutate(pool = factor(pool,
                             levels = c("soil", "detritus", "plant",
                                        "detritivore", "herbivore", "predator"))) %>%
        arrange(pool, time) %>%
        group_by(pool) %>%
        summarize(max_gain_V = max(Vg - Vg[1]),
                  min_gain_V = min(Vg - Vg[1]),
                  max_loss_V = max(Vl - Vl[1]),
                  min_loss_V = min(Vl - Vl[1]),
                  max_gain_H = max(Hg - Hg[1]),
                  min_gain_H = min(Hg - Hg[1]),
                  max_loss_H = max(Hl - Hl[1]),
                  min_loss_H = min(Hl - Hl[1]),
                  cum_pos_loss_V = sum(Vl[Vl > Vl[1]] - Vl[1]),
                  cum_neg_loss_V = sum(Vl[Vl < Vl[1]] - Vl[1]),
                  cum_pos_loss_H = sum(Hl[Hl > Hl[1]] - Hl[1]),
                  cum_neg_loss_H = sum(Hl[Hl < Hl[1]] - Hl[1]),
                  cum_gain_V = sum(Vg[Vg > Vg[1]] - Vg[1]),
                  cum_gain_H = sum(Hg[Hg > Hg[1]] - Hg[1]),

                  cum_MfR = sum(MfR),
                  cum_MfD = sum(MfD)) %>%
        ungroup() %>%
        mutate(w = .w, b = .b, f = .other_pars$f, area = b * w) %>%
        select(w, b, f, area, everything())
    return(fw)
}


# Equilibrium pool size changes
ep_df <- crossing(aDV = c(par_estimates$aDV[1], par_estimates$aPH[1]),
                  aPH = c(par_estimates$aPH[1], par_estimates$aDV[1])) %>%
    mutate(pools = map2(aDV, aPH, ~ equil_pools(aDV = .x, aPH = .y)))




# ------------------------
# Combinations of parameter values
# ------------------------

par_combs <- expand.grid(w = seq(10, 25, length.out = 25),
                         b = seq(0.1, 40, length.out = 25),
                         f = seq(8e-3, 8e-2, length.out = 25),
                         mM = par_estimates$mM[1] * c(0.5, 1, 2),
                         hM = par_estimates$hM[1] * c(0.5, 1, 2),
                         # Plant / herbivore uptake rates:
                         aDV = c(par_estimates$aDV[1], par_estimates$aPH[1]),
                         aPH = c(par_estimates$aPH[1], par_estimates$aDV[1])) %>%
    split(row(.)[,1])


# ------------------------
# Run and summarize simulations
# ------------------------

pulse_df <- pbmclapply(par_combs, one_combo, mc.cores = n_cores)
pulse_df <- bind_rows(pulse_df)


write_csv(pulse_df, "data-raw/pulse_data.csv")

