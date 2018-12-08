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
                      c(N = iN - aNP*N*P*(1-P/kP) + (1-lD)*mD*D - mN*N,
                        D = (1-lP)*mP*P + (1-lV)*mV*V + (1-lH)*mH*H +
                            (1-lR)*mR*R + (1 - lM)*mM*M - aDV*D*V*(1-V/kV) -
                            mD*D,
                        P = aNP*N*P*(1-P/kP) - aPH*P*H*(1-H/kH) - mP*P,
                        V = aDV*D*V*(1-V/kV) - (aR*V*R)*(1-R/kR) - mV*V,
                        H = aPH*P*H*(1-H/kH) - (aR*H*R)*(1-R/kR) - mH*H,
                        R = (aR*V*R + aR*H*R + aR*M*R)*(1-R/kR) - mR*R,
                        M = iM - mM*M - (aR*M*R)*(1-R/kR))
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
#' @param r Period of the pulse expresed in units of \code{1.5 * w} (so the pulses
#'     don't overlap). Also controls the smoothness, along with \code{a}.
#' @param w Width of the midge pulse.
#' @param d Mid point of the first pulse in units of \code{w}.
#' @param tstep Step size in units of time. Defaults to \code{1}.
#' @param V Boolean for whether to include the V pool.
#' @param R Boolean for whether to include the R pool.
#'     This pool isn't allowed when no V or H pool are present.
#' @param H Boolean for whether to include the H pool.
#' @param pool_starts Named list of initial values for each pool.
#'     Names can include "N0", "D0", "P0", "V0", "H0", "R0", and "M0".
#'     Any pools not included here will start at their equilibrium value, except
#'     for midges that start with zero by default.
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
#'
#' @export
#'
#'
food_web <- function(tmax, a, b, r, w, d, tstep = 1, V = TRUE, R = TRUE, H = TRUE,
                     pool_starts = NULL) {

    if (!inherits(V, "logical") || !inherits(R, "logical") || !inherits(H, "logical") ||
        length(V) != 1 || length(R) != 1 || length(H) != 1) {
        stop("\nV, R, and H must all be length-1 logicals (TRUE or FALSE).")
    }
    if (!H && !V && R) {
        stop("\nYou can't have predators if you have no detritivores or herbivores, ",
             "you silly goose.")
    }
    if (!V || !R || !H) {
        stop("\nFood webs other than the full one aren't yet programmed. (Blame Joe.)")
    }

    pars <- par_estimates %>%
        filter(V == V, R == R, H == H) %>%
        dplyr::select(-V, -R, -H) %>%
        {attr(., "spec") <- NULL; .} %>%
        as.list()
    pars$midges <- function(t_) midge_pulse(t_, a, b, r, w, d)

    pool_names <- c("N", "D", "P", "V", "H", "R", "M")
    # Default values:
    init <- c(unlist(pars[paste0(pool_names[pool_names != "M"], "eq")]), M0 = 0)
    names(init) <- c(paste0(pool_names[pool_names != "M"], "0"), "M0")
    # Replace those that are provided:
    if (!is.null(pool_starts)) {
        if (!inherits(pool_starts, "list")) {
            stop("\nIf not NULL, pool_starts must be a list.")
        }
        if (is.null(names(pool_starts)) || any(is.na(names(pool_starts)))) {
            stop("\nIn pool_starts, all items must be named.");
        }
        if (any(!names(pool_starts) %in% paste0(pool_names, "0"))) {
            stop("\nThe following name(s) in pool_starts don't match with a pool ",
                 "starting-value argument name: ",
                 names(pool_starts)[!names(pool_starts) %in% paste0(pool_names, "0")])
        }
        init[names(pool_starts)] <- unlist(pool_starts)
    }
    names(init) <- c(pool_names[pool_names != "M"], "M")
    if (!V) init$V0 <- 0
    if (!R) init$R0 <- 0
    if (!H) init$H0 <- 0

    solved_ode = ode(init, seq(0, tmax, tstep), diff_eq, pars)

    solved_ode <- solved_ode %>%
        as.data.frame() %>%  # prevents "matrix as column is not supported" error
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

    ;

}
