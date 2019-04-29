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
                      c(N = iN - aNP*N*P/(1 + aNP*hN*N) + (1 - lD)*mD*D - mN*N,
                        D = (1 - lP)*(mP0 + mP*P)*P + (1 - lV)*(mV0 + mV*V)*V +
                            (1 - lH)*(mH0 + mH*H)*H + (1 - lR)*(mR0 + mR*R)*R +
                            (1 - lM)*mM*M - aDV*D*V/(1 + aDV*hD*D) - mD*D,
                        P = aNP*N*P/(1 + aNP*hN*N) - aPH*P*H/(1 + aPH*hP*P) - (mP0 + mP*P)*P,
                        V = aDV*D*V/(1 + aDV*hD*D) - (aR*V*R)/(1 + aR*hVHM*(V + H) + (aR * f)*hVHM*M) - (mV0 + mV*V)*V,
                        H = aPH*P*H/(1 + aPH*hP*P) - (aR*H*R)/(1 + aR*hVHM*(V + H) + (aR * f)*hVHM*M) - (mH0 + mH*H)*H,
                        R = (aR*V*R + aR*H*R + (aR * f)*M*R)/(1 + aR*hVHM*(V + H) + (aR * f)*hVHM*M) - (mR0 + mR*R)*R,
                        M = iM - mM*M - ((aR * f)*M*R)/(1 + aR*hVHM*(V + H) + (aR * f)*hVHM*M))

                  }))

    return(list(output))
}


#' Run a food web.
#'
#' For more info on the midge-pulse parameters (`b`, `s`, `w`),
#' see `vignette("smooth_pulse", "mtf")`.
#'
#'
#' @param tmax Duration over which to run the model.
#' @param b Maximum value of the midge pulse.
#' @param s Start time for midge pulse
#' @param w Width of the midge pulse.
#' @param tstep Step size in units of time. Defaults to \code{1}.
#' @param pool_starts Named list of initial values for each pool.
#'     Names can include "N0", "D0", "P0", "V0", "H0", "R0", and "M0"
#'     (for nitrogen, detritus, plant, detritivore, herbivore, predator, and midge
#'     pools, respectively).
#'     Any pools not included here will start at their equilibrium value, except
#'     for midges that start at zero by default.
#' @param ep_obj Output from the \code{\link{equil_pools}} function that specifies
#'     equilibrium pool sizes when you change variable values from defaults.
#'     This is used when you pass something to `other_pars` that changes equilibrium
#'     pool size(s).
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
#' @importFrom purrr set_names
#'
#' @return A data frame with the following columns: `time`, `pool`, `N`.
#'
#'
#'
#' @export
#'
#' @examples
#' # What happens to the food web bc of this midge pulse?
#' web_output <- food_web(tmax = 1000, b = 400, s = 10, w = 40)
#' web_output
#'
#'
food_web <- function(tmax, b, s, w, tstep = 1,
                     pool_starts = NULL,
                     ep_obj = NULL,
                     other_pars = list()) {

    .iN <- 10
    if (!is.null(other_pars$iN)) .iN <- other_pars$iN

    pars <- par_estimates %>%
        filter(V == 1, R == 1, H == 1, iN == .iN) %>%
        dplyr::select(-V, -R, -H)
    if (nrow(pars) == 0) {
        iN_vals <- par_estimates %>%
            filter(V == 1, R == 1, H == 1) %>%
            .[["iN"]] %>%
            unique() %>%
            paste(collapse = ", ")
        stop("\nYou're requesting an `iN` value that we don't have parameter values ",
             "for. The possibilities are as follows: ", iN_vals)
    }
    pars <- as.list(pars)
    pars$midges <- function(t_) midge_pulse(t_, b, s, w)
    if (length(other_pars) > 0) {
        if (!all(names(other_pars) %in% names(pars))) {
            bad_pn <- names(other_pars)[!names(other_pars) %in% names(pars)]
            stop("\nThe following name(s) in ... args don't match with a parameter ",
                 "for the simulations: ",
                 paste(sprintf('"%s"', bad_pn), collapse = ", "))
        }
        for (p in names(other_pars)) pars[[p]] <- other_pars[[p]]
    }

    # Set equilibrium values from `ep_obj`
    if (!is.null(ep_obj)) {
        stopifnot(inherits(ep_obj, "equil_pools"))
        ep_vals <- ep_obj$vals %>%
            .[["eq"]] %>%
            set_names(nm = paste0(c("N", "D", "P", "V", "H", "R"), "eq"))
        for (p in names(ep_vals)) pars[[p]] <- ep_vals[[p]]
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

    time <- seq(0, tmax, tstep)
    if (time[length(time)] < tmax) time <- c(time, tmax)
    solved_ode <- ode(init, time, diff_eq, pars)

    solved_ode <- solved_ode %>%
        as.data.frame() %>%  # <-- prevents "matrix as column is not supported" error
        as_tibble() %>%
        rename(
            soil = N,
            detritus = D,
            plant = P,
            detritivore = V,
            herbivore = H,
            predator = R,
            midge = M
        ) %>%
        gather('pool', 'N', -time) %>%
        mutate(pool = factor(pool, levels = c("soil", "detritus", "plant",
                                              "detritivore", "herbivore",
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
        geom_hline(data = web_df %>% filter(time == min(time)),
                   aes(yintercept=N), color="firebrick") +
        geom_line(size = 1) +
        geom_point(aes(time, minb), shape="") +
        theme_classic()
}
