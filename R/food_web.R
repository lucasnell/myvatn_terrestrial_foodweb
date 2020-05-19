fI <- function(I, P, D, pars) {
    pars$M <- 0
    with(pars, iI - aIP*I*P/(1 + aIP*hI*I) + (1 - lD)*muD*D - muI*I)
}
fD <- function(D, P, V, H, X, pars) {
    pars$M <- 0
    with(pars, (1 - lP)*(muP + mP*P)*P + (1 - lV)*(muV + mV*V)*V + (1 - lH)*(muH + mH*H)*H + (1 - lX)*(muX + mX*X)*X + (1 - lM)*muM*M - aDV*D*V/(1 + aDV*hD*D) - muD*D)
}
fP <- function(P, I, H, pars) {
    pars$M <- 0
    with(pars, aIP*I*P/(1 + aIP*hI*I) - aPH*P*H/(1 + aPH*hP*P) - (muP + mP*P)*P)
}
fV <- function(V, D, X, H, pars) {
    pars$M <- 0
    with(pars, aDV*D*V/(1 + aDV*hD*D) - (aX*V*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M) - (muV + mV*V)*V)
}
fH <- function(H, P, X, V, pars) {
    pars$M <- 0
    with(pars, aPH*P*H/(1 + aPH*hP*P) - (aX*H*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M) - (muH + mH*H)*H)
}
fX <- function(X, V, H, pars) {
    pars$M <- 0
    with(pars, (aX*V*X + aX*H*X + (aX * q)*M*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M) - (muX + mX*X)*X)
}



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
                      c(I = iI - aIP*I*P/(1 + aIP*hI*I) + (1 - lD)*muD*D - muI*I,
                        D = (1 - lP)*(muP + mP*P)*P + (1 - lV)*(muV + mV*V)*V +
                            (1 - lH)*(muH + mH*H)*H + (1 - lX)*(muX + mX*X)*X +
                            (1 - lM)*muM*M - aDV*D*V/(1 + aDV*hD*D) - muD*D,
                        P = aIP*I*P/(1 + aIP*hI*I) - aPH*P*H/(1 + aPH*hP*P) - (muP + mP*P)*P,
                        V = aDV*D*V/(1 + aDV*hD*D) - (aX*V*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M) - (muV + mV*V)*V,
                        H = aPH*P*H/(1 + aPH*hP*P) - (aX*H*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M) - (muH + mH*H)*H,
                        X = (aX*V*X + aX*H*X + (aX * q)*M*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M) - (muX + mX*X)*X,
                        M = iM - muM*M - ((aX * q)*M*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M))
                  }))

    return(list(output))
}


#' Differential equations for food web.
#'
#' In this version, midges DON'T go to predators!
#'
#' @noRd
#'
diff_eq__no_M_to_X <- function(t, y, pars) {

    iM = pars$midges(t)

    output <-
        with(as.list(pars),
             with(as.list(y),
                  {
                      c(I = iI - aIP*I*P/(1 + aIP*hI*I) + (1 - lD)*muD*D - muI*I,
                        D = (1 - lP)*(muP + mP*P)*P + (1 - lV)*(muV + mV*V)*V +
                            (1 - lH)*(muH + mH*H)*H + (1 - lX)*(muX + mX*X)*X +
                            (1 - lM)*muM*M - aDV*D*V/(1 + aDV*hD*D) - muD*D,
                        P = aIP*I*P/(1 + aIP*hI*I) - aPH*P*H/(1 + aPH*hP*P) - (muP + mP*P)*P,
                        V = aDV*D*V/(1 + aDV*hD*D) - (aX*V*X)/(1 + aX*hX*(V + H)) - (muV + mV*V)*V,
                        H = aPH*P*H/(1 + aPH*hP*P) - (aX*H*X)/(1 + aX*hX*(V + H)) - (muH + mH*H)*H,
                        X = (aX*V*X + aX*H*X)/(1 + aX*hX*(V + H)) - (muX + mX*X)*X,
                        M = iM - muM*M)
                  }))

    return(list(output))

}



#' Differential equations for food web.
#'
#' In this version, midges DON'T go to detritus!
#'
#' @noRd
#'
diff_eq__no_M_to_D <- function(t, y, pars) {

    iM = pars$midges(t)

    output <-
        with(as.list(pars),
             with(as.list(y),
                  {
                      c(I = iI - aIP*I*P/(1 + aIP*hI*I) + (1 - lD)*muD*D - muI*I,
                        D = (1 - lP)*(muP + mP*P)*P + (1 - lV)*(muV + mV*V)*V +
                            (1 - lH)*(muH + mH*H)*H + (1 - lX)*(muX + mX*X)*X -
                            aDV*D*V/(1 + aDV*hD*D) - muD*D,
                        P = aIP*I*P/(1 + aIP*hI*I) - aPH*P*H/(1 + aPH*hP*P) - (muP + mP*P)*P,
                        V = aDV*D*V/(1 + aDV*hD*D) - (aX*V*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M) - (muV + mV*V)*V,
                        H = aPH*P*H/(1 + aPH*hP*P) - (aX*H*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M) - (muH + mH*H)*H,
                        X = (aX*V*X + aX*H*X + (aX * q)*M*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M) - (muX + mX*X)*X,
                        M = iM - muM*M - ((aX * q)*M*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M))
                  }))

    return(list(output))

}




#' Run a food web.
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
#' @param midges_not_to Character specifying whether midges should NOT go to detritus
#'     or NOT to predators. If `"none"`, `NULL`, or `NA`, midges go to both.
#'     Takes as input the following: `"none"`, `NULL`, `NA`, `"X"`, `"predator"`,
#'     `"D"`, or `"detritus"`. Defaults to `NULL`.
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
                     other_pars = list(),
                     midges_not_to = NULL) {

    .iI <- 10
    if (!is.null(other_pars$iI)) .iI <- other_pars$iI

    if (!is.null(midges_not_to) && !is.na(midges_not_to)) {
        midges_not_to <- match.arg(midges_not_to,
                                   c("none", "X", "predator", "D", "detritus"))
    }

    pars <- par_estimates %>%
        filter(V == 1, X == 1, H == 1, iI == .iI) %>%
        dplyr::select(-V, -X, -H)
    if (nrow(pars) == 0) {
        iI_vals <- par_estimates %>%
            filter(V == 1, X == 1, H == 1) %>%
            .[["iI"]] %>%
            unique() %>%
            paste(collapse = ", ")
        stop("\nYou're requesting an `iI` value that we don't have parameter values ",
             "for. The possibilities are as follows: ", iI_vals)
    }
    pars <- as.list(pars)
    pars$midges <- function(t_) {
        i_M_t = ifelse(t_ > s & t_ <= (s + w), b, 0)
        return(i_M_t)
    }
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
            set_names(nm = paste0(c("I", "D", "P", "V", "H", "X"), "eq"))
        for (p in names(ep_vals)) pars[[p]] <- ep_vals[[p]]
    }

    pool_names <- c("I", "D", "P", "V", "H", "X", "M")
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

    .diff_eq <- diff_eq
    if (isTRUE(midges_not_to %in% c("X", "predator"))) .diff_eq <- diff_eq__no_M_to_X
    if (isTRUE(midges_not_to %in% c("D", "detritus"))) .diff_eq <- diff_eq__no_M_to_D
    solved_ode <- ode(init, time, .diff_eq, pars)

    solved_ode <- solved_ode %>%
        as.data.frame() %>%  # <-- prevents "matrix as column is not supported" error
        as_tibble() %>%
        rename(
            soil = I,
            detritus = D,
            plant = P,
            detritivore = V,
            herbivore = H,
            predator = X,
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


