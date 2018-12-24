#' Differential equations for food web.
#'
#'
#' @noRd
#'
diff_eq <- function(t, y, pars) {

    iM = pars$midges(t)

    output <-
        with(as.list(pars),
             with(as.list(y),
                  {
                      c(N = iN - aNP*N*P/(hN + N) + (1 - lD)*mD*D - mN*N,
                        D = (1 - lP)*(mP0 + mP*P)*P + (1 - lV)*(mV0 + mV*V)*V +
                            (1 - lH)*(mH0 + mH*H)*H + (1 - lR)*(mR0 + mR*R)*R +
                            (1 - lM)*mM*M - aDV*D*V/(hD + D) - mD*D,
                        P = aNP*N*P/(hN + N) - aPH*P*H/(hP + P) - (mP0 + mP*P)*P,
                        V = aDV*D*V/(hD + D) - (aR*V*R)/(hVHM + V + H + M) - (mV0 + mV*V)*V,
                        H = aPH*P*H/(hP + P) - (aR*H*R)/(hVHM + V + H + M) - (mH0 + mH*H)*H,
                        R = (aR*V*R + aR*H*R + aM*M*R)/(hVHM + V + H + M) - (mR0 + mR*R)*R,
                        M = iM - mM*M - (aM*M*R)/(hVHM + V + H + M))
                  }))

    return(list(output))
}


#' Run a food web.
#'
#' For more info on the midge-pulse parameters (`a`, `b`, `r`, `w`, and `d`),
#' see `vignette("smooth_pulse", "mtf")`.
#'
#'
#' @param tmax Duration over which to run the model.
#' @param a Controls the smoothness of the pulse, along with \code{r}. If \code{a} is
#'     sufficiently high, the pulse will always be rectangular.
#' @param b Maximum value of the midge pulse.
#' @param r Period of the pulse expressed in units of \code{1.5 * w} (so the pulses
#'     don't overlap). Also controls the smoothness, along with \code{a}.
#' @param w Width of the midge pulse.
#' @param d Mid point of the first pulse in units of \code{w}.
#' @param tstep Step size in units of time. Defaults to \code{1}.
#' @param .V Boolean for whether to include the V pool. Defaults to `TRUE`.
#' @param .R Boolean for whether to include the R pool. Defaults to `TRUE`.
#' @param .H Boolean for whether to include the H pool. Defaults to `TRUE`.
#' @param .iN Input rate from pool N. Options for this include `200`, `1000`, or `1800`.
#'     Defaults to `1000`.
#' @param pool_starts Named list of initial values for each pool.
#'     Names can include "N0", "D0", "P0", "V0", "H0", "R0", and "M0"
#'     (for nitrogen, detritus, plant, detrivore, herbivore, predator, and midge
#'     pools, respectively).
#'     Any pools not included here will start at their equilibrium value, except
#'     for midges that start at zero by default.
#' @param other_pars A named list of other parameter values to use instead of defaults.
#'     Possible parameter names include all column names in `par_estimates`
#'     except the first three (see `?par_estimates` for a description of column names).
#'
#'
#'
#' @importFrom deSolve ode
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr as_data_frame
#' @importFrom dplyr rename
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom tidyr gather
#'
#' @return A data frame with the following columns: `time`, `pool`, `N`.
#'
#'
#' @usage food_web(tmax, a, b, r, w, d,
#'          tstep = 1,
#'          .V = TRUE, .R = TRUE, .H = TRUE,
#'          .iN = 1000,
#'          pool_starts = NULL)
#'
#' @export
#'
#' @examples
#' # Test a midge pulse:
#' midges <- test_midges(1000, a = 1e9, b = 10, r = 400, w = 40, d = 1)
#' plot(midges, type = 'l', ylab = "iM")
#'
#' # What happens to the food web bc of this midge pulse?
#' web_output <- food_web(tmax = 1000, a = 1e9, b = 0, r = 400, w = 40, d = 1)
#' web_output
#'
#' # If you want to adjust starting value for predators:
#' web_output <- food_web(tmax = 1000, a = 1e9, b = 0, r = 400, w = 40, d = 1,
#'                        pool_starts = list(R0 = 10))
#' web_output
#'
#' # If you want to remove predators from the food web entirely:
#' web_output <- food_web(tmax = 1000, a = 1e9, b = 0, r = 400, w = 40, d = 1,
#'                        .R = FALSE)
#' web_output
#'
#' # Plotting output:
#' plot_web(web_output)
#'
food_web <- function(tmax, a, b, r, w, d, tstep = 1,
                     .V = TRUE, .R = TRUE, .H = TRUE,
                     .iN = 1000,
                     pool_starts = NULL,
                     other_pars = NULL) {

    if (!inherits(.V, "logical") || !inherits(.R, "logical") || !inherits(.H, "logical") ||
        length(.V) != 1 || length(.R) != 1 || length(.H) != 1) {
        stop("\n.V, .R, and .H must all be length-1 logicals (TRUE or FALSE).")
    }
    if (!.iN %in% par_estimates$iN) {
        stop("\nYour options for .iN are ",
             paste(unique(par_estimates$iN), collapse = ", "))
    }

    pars <- par_estimates %>%
        filter(V == ifelse(.V,1,0), R == ifelse(.R,1,0), H == ifelse(.H,1,0),
               iN == .iN) %>%
        dplyr::select(-V, -R, -H)
    if (nrow(pars) == 0) {
        stop("\nYou're using a combination of .V, .R, .H, and .iN that we don't have ",
             "parameter values for.")
    }
    pars <- as.list(pars)
    pars$midges <- function(t_) midge_pulse(t_, a, b, r, w, d)
    if (!is.null(other_pars)) {
        if (!inherits(other_pars, "list") || is.null(names(other_pars))) {
            stop("other_pars must be NULL or a named list.")
        }
        if (!all(names(other_pars) %in% colnames(par_estimates)[-1:-3])) {
            bad_pars <- names(other_pars)[!names(other_pars) %in%
                                              colnames(par_estimates)[-1:-3]]
            stop("The following name(s) in other_pars aren't present in ",
                 "`colnames(par_estimates)[-1:-3]`: ",
                 paste(sprintf('"%s"', bad_pars), collapse = ", "), ".")
        }
        pars[names(other_pars)] <- other_pars
    }

    pool_names <- c("N", "D", "P", "V", "H", "R", "M")
    # Default values:
    init <- c(unlist(pars[paste0(pool_names[pool_names != "M"], "eq")]), M0 = 0)
    names(init) <- c(paste0(pool_names[pool_names != "M"], "0"), "M0")
    # Replace those that are provided:
    if (!is.null(pool_starts)) {
        if (!inherits(pool_starts, "list") || is.null(names(pool_starts))) {
            stop("\nIf not NULL, pool_starts must be a named list.")
        }
        if (any(!names(pool_starts) %in% paste0(pool_names, "0"))) {
            bad_ps <- names(pool_starts)[!names(pool_starts) %in% paste0(pool_names, "0")]
            stop("\nThe following name(s) in pool_starts don't match with a pool ",
                 "starting-value argument name: ",
                 paste(sprintf('"%s"', bad_ps), collapse = ", "))
        }
        init[names(pool_starts)] <- unlist(pool_starts)
    }
    names(init) <- c(pool_names[pool_names != "M"], "M")
    if (!.V) init[["V"]] <- 0
    if (!.R) init[["R"]] <- 0
    if (!.H) init[["H"]] <- 0

    time <- seq(0, tmax, tstep)
    if (time[length(time)] < tmax) time <- c(time, tmax)
    solved_ode <- ode(init, time, diff_eq, pars)

    solved_ode <- solved_ode %>%
        as.data.frame() %>%  # <-- prevents "matrix as column is not supported" error
        as_data_frame() %>%
        rename(
            nitrogen = N,
            detritus = D,
            plant = P,
            detrivore = V,
            herbivore = H,
            predator = R,
            midge = M
        ) %>%
        gather('pool', 'N', -time) %>%
        mutate(pool = factor(pool, levels = c("nitrogen", "detritus", "plant",
                                              "detrivore", "herbivore",
                                              "predator", "midge")),
               time = as.integer(time)) %>%
        arrange(pool, time)

    return(solved_ode)

}



#' Plot output from the `food_web` function.
#'
#' @param web_df Dataframe output from the `food_web` function.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @import ggplot2
#'
#' @export
#'
#'
plot_web <- function(web_df) {
    web_df %>%
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
}
