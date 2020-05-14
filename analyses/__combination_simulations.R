
# Run the simulations analyzed in long-term_dynamics.R


if (!require(pbmcapply)) install.packages("pbmcapply")


# load packages
suppressPackageStartupMessages({
    library(mtf)
    library(tidyverse)
    library(pbmcapply)
})

# Number of threads to use:
n_threads <- as.integer(detectCores() * 0.9)

# ------------------------
# Summary functions
# ------------------------

parlist <- par_estimates %>%
        filter(V==1, H==1, X==1, iI == 10) %>%
        as.list()
V_gain <- function(V, D, aDV) {
    hD <- parlist[["hD"]]
    (aDV*D*V/(1 + aDV*hD*D)) / V
}
V_loss <- function(V, X, H, M, q, hM) {
    aX <- parlist[["aX"]]
    hX <- parlist[["hX"]]
    ((aX*V*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M)) / V
}
H_gain <- function(P, H, aPH) {
    hP <- parlist[["hP"]]
    (aPH*P*H/(1 + aPH*hP*P)) / H
}
H_loss <- function(H, X, V, M, q, hM) {
    aX <- parlist[["aX"]]
    hX <- parlist[["hX"]]
    ((aX*H*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M)) / H
}

M_flux_X <- function(H, X, V, M, q, hM) {
    aX <- parlist[["aX"]]
    hX <- parlist[["hX"]]
    ((aX * q)*M*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M) / M
}
M_flux_D <- function(M, muM) {
    muM*M
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

    .q <- ifelse(!is.null(.other_pars$q), .other_pars$q, parlist$q)
    .muM <- ifelse(!is.null(.other_pars$muM), .other_pars$muM, parlist$muM)
    .hM <- ifelse(!is.null(.other_pars$hM), .other_pars$hM, parlist$hM)
    .aDV <- ifelse(!is.null(.other_pars$aDV), .other_pars$aDV, parlist$aDV)
    .aPH <- ifelse(!is.null(.other_pars$aPH), .other_pars$aPH, parlist$aPH)


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

    fw <- food_web(tmax = 500, s = 10, b = .b, w = .w, ep_obj = .ep_obj,
                   other_pars = .other_pars)

    fw <- fw %>%
        spread(pool, N) %>%
        mutate(Vg = V_gain(detritivore, detritus, .other_pars$aDV),
               Vl = V_loss(detritivore, predator, herbivore, midge, .other_pars$q,
                           .other_pars$hM),
               Hg = H_gain(plant, herbivore, .other_pars$aPH),
               Hl = H_loss(herbivore, predator, detritivore, midge, .other_pars$q,
                           .other_pars$hM),
               MfX = M_flux_X(herbivore, predator, detritivore, midge, .other_pars$q,
                        .other_pars$hM),
               MfD = M_flux_D(midge, .other_pars$muM)) %>%
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

                  cum_MfX = sum(MfX),
                  cum_MfD = sum(MfD)) %>%
        ungroup() %>%
        mutate(w = .w,
               b = .b,
               q = .q,
               muM = .muM,
               hM = .hM,
               aDV = .aDV,
               aPH = .aPH) %>%
        select(w, b, q, muM, hM, aDV, aPH, everything())
    return(fw)
}


# Equilibrium pool size changes
ep_df <- crossing(aDV = c(par_estimates$aDV[1], par_estimates$aPH[1]),
                  aPH = c(par_estimates$aPH[1], par_estimates$aDV[1])) %>%
    mutate(pools = map2(aDV, aPH, ~ equil_pools(aDV = .x, aPH = .y)))










# ------------------------
# Combinations of parameter values
# ------------------------

par_combs <- expand.grid(w = 20,
                         b = seq(0.1, 50, length.out = 25),
                         q = c(0.8, 4.4, 8),
                         muM = par_estimates$muM[1] * c(0.5, 1, 2),
                         hM = par_estimates$hM[1] * c(0.5, 1, 2),
                         # Plant / herbivore uptake rates:
                         aDV = c(par_estimates$aDV[1], par_estimates$aPH[1]),
                         aPH = c(par_estimates$aPH[1], par_estimates$aDV[1])) %>%
    # This pairing makes detritivores go to zero:
    filter(!(aDV == par_estimates$aPH[1] & aPH == par_estimates$aDV[1])) %>%
    split(row(.)[,1])


# ------------------------
# Run and summarize simulations
# ------------------------

pulse_df <- pbmclapply(par_combs, one_combo, mc.cores = n_threads)
pulse_df <- bind_rows(pulse_df)

write_csv(pulse_df, "~/Box Sync/Iceland Food Web Model/Results/sim_combinations.csv.gz")

