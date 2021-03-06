
#' Calculate equilibrium pool sizes after changing parameters.
#'
#' This function only works for when all pools are present.
#'
#' @inheritParams food_web
#' @param ... Other parameters to calculate equilibrium pool sizes for.
#'     All parameters in \code{\link{par_estimates}} other than `V`, `H`, and `X`, plus
#'     those ending in `eq` are permitted.
#'
#' @importFrom deSolve ode
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr as_tibble
#' @importFrom dplyr rename
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom tidyr gather
#'
#' @export
#'
equil_pools <- function(tmax = 1000, tstep = 1,
                        midges_not_to = NULL,
                        ...) {


    .call <- match.call()

    if (!is.null(midges_not_to) && !is.na(midges_not_to)) {
        midges_not_to <- match.arg(midges_not_to,
                                   c("none", "X", "predator", "D", "detritus"))
    }

    pars <- par_estimates %>%
        filter(V == 1, X == 1, H == 1, iI == iI[2]) %>%
        dplyr::select(-V, -X, -H) %>%
        as.list()

    pars$midges <- function(t_) return(0)

    other_pars <- list(...)
    if (length(other_pars) > 0) {
        if (!all(names(other_pars) %in% names(pars))) {
            bad_pn <- names(other_pars)[!names(other_pars) %in% names(pars)]
            stop("\nThe following name(s) in ... args don't match with a parameter ",
                 "for the simulations: ",
                 paste(sprintf('"%s"', bad_pn), collapse = ", "))
        }
        for (p in names(other_pars)) pars[[p]] <- other_pars[[p]]
    }

    # Starting pool sizes (no reason for user to adjust these directly):
    pool_names <- c("I", "D", "P", "V", "H", "X", "M")
    init <- c(unlist(pars[paste0(pool_names[pool_names != "M"], "eq")]), M0 = 0)
    names(init) <- pool_names


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
        filter(pool != "midge") %>%
        mutate(pool = factor(pool, levels = c("soil", "detritus", "plant",
                                              "detritivore", "herbivore",
                                              "predator")),
               time = as.integer(time)) %>%
        arrange(pool, time)

    eq_pools <- solved_ode %>%
        filter(time == tmax) %>%
        group_by(pool) %>%
        summarize(eq = N) %>%
        ungroup()

    out <- structure(list(vals = eq_pools, ts = solved_ode, call = .call),
                     class = "equil_pools")

    return(out)

}



#' Plot time series of `equil_pools` object.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 scale_color_brewer
#'
#' @export
#' @noRd
#'
plot.equil_pools <- function(x, ...) {

    x$ts %>%
        filter(pool != "midge") %>%
        ggplot(aes(time, N)) +
        geom_line(aes(color = pool)) +
        facet_wrap(~ pool, scales = "free_y") +
        scale_color_brewer(palette = "Dark2", guide = FALSE)
}



#' Print `equil_pools` object.
#'
#' @export
#' @noRd
#'
print.equil_pools <- function(x, ...) {

    cat("<New equilibrium values>\n\n")
    print(x$vals)
    invisible(NULL)
}
