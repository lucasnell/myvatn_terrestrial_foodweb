
#' Foodweb class that implements either of the two versions of the model (A and B).
#'
#' Upon creating a \code{web} object, it solves for selected unknown parameters, given
#' the 'known' parameters and equilibria.
#' It also includes an internal class function \code{ode_solve} that solves the ODEs.
#'
#'
#' @usage
#' \preformatted{w <- web$new(...)
#'
#' w$eq_solve(solve_pars = c("mP", "mD", "aNP", "aDV", "aPH", "aR"),
#'            initial_vals = rep(0.1, 6))
#' w$values()
#' w$ode_solve(tmax, a, b, r, w, d, tstep = 1)
#' w$test_midges(tmax, a, b, r, w, d, tstep = 1)
#' print(w)
#' }
#'
#'
#'
#' @param ... Custom settings for any of the slots. See below for possible arguments.
#' @param solve_pars A character vector of length 6 of the parameters that need to be
#'     solved for. Defaults to \code{c("mP", "mD", "aNP", "aDV", "aPH", "aR")}.
#' @param initial_vals Numeric vector of length 6 representing the initial guesses of
#'     the unknown parameters. For more info, see \code{\link[rootSolve]{multiroot}}.
#' @param tmax Duration over which to run the model.
#' @param a Controls the smoothness of the pulse, along with \code{r}. If \code{a} is
#'     sufficiently high, the pulse will always be rectangular (for more info,
#'     see `vignette("smooth_pulse", "mtf")`).
#' @param b Maximum value of the midge pulse (for more info, see
#'     `vignette("smooth_pulse", "mtf")`).
#' @param r Period of the pulse expresed in units of \eqn{1.5 \times w} (so the pulses
#'     don't overlap). Also controls the smoothness, along with \code{a} (for more
#'     info, see `vignette("smooth_pulse", "mtf")`).
#' @param w Width of the midge pulse (for more info, see
#'     `vignette("smooth_pulse", "mtf")`).
#' @param d Mid point of the first pulse in units of \eqn{w} (for more info, see
#'     `vignette("smooth_pulse", "mtf")`).
#' @param tstep Step size in units of time. Defaults to \code{1}.
#'
#'
#'
#' @section Methods:
#'
#' \describe{
#'     \item{\code{$values()}}{
#'         Return a list of all parameter values.
#'     }
#'     \item{\code{$eq_solve(solve_pars, initial_vals)}}{
#'         Use estimated equilibrium biomass and selected values for certain
#'         parameters to solve for unknown values. "Known" parameters must be selected
#'         with care, to ensure that an equilibrium solution can actually be reached.
#'     }
#'     \item{\code{$ode_solve(tmax, a, b, r, w, d, tstep)}}{
#'         Solve the ODE and output timeseries of nitrogen content for each pool.
#'     }
#'     \item{\code{$test_midges(tmax, a, b, r, w, d, tstep)}}{
#'         Test out a particular midge pulse scenario. It returns a vector of
#'         inputs of midges (in units of N) to the midge pool through time.
#'     }
#' }
#'
#'
#'
#' @section Slots:
#'
#' @slot Neq  Desired equilibrium value for N. Initially set to \code{34300}.
#' @slot Deq  Desired equilibrium value for D. Initially set to \code{308700}.
#' @slot Peq  Desired equilibrium value for P. Initially set to \code{4300}.
#' @slot Veq  Desired equilibrium value for V. Initially set to \code{81}.
#' @slot Heq  Desired equilibrium value for H. Initially set to \code{24}.
#' @slot Req  Desired equilibrium value for R. Initially set to \code{13}.
#' @slot Meq  Desired equilibrium value for M. Initially set to \code{0}.
#' @slot N0 Initial state for N. Initially set to \code{NULL}.
#' @slot D0 Initial state for D. Initially set to \code{NULL}.
#' @slot P0 Initial state for P. Initially set to \code{NULL}.
#' @slot V0 Initial state for V. Initially set to \code{NULL}.
#' @slot H0 Initial state for H. Initially set to \code{NULL}.
#' @slot R0 Initial state for R. Initially set to \code{NULL}.
#' @slot M0 Initial state for M. Initially set to \code{NULL}.
#' @slot iN Input to N. Initially set to \code{1000}.
#' @slot lD Loss rates systems for D. Initially set to \code{0.1}.
#' @slot lP Loss rates systems for P. Initially set to \code{0.1}.
#' @slot lV Loss rates systems for V. Initially set to \code{0.1}.
#' @slot lH Loss rates systems for H. Initially set to \code{0.1}.
#' @slot lR Loss rates systems for R. Initially set to \code{0.1}.
#' @slot lM Loss rates systems for M. Initially set to \code{0.1}.
#' @slot mN Loss rates from pool N (returned to either N or D).
#'     Initially set to \code{0.002}.
#' @slot mP Loss rates from pool P (returned to either N or D).
#'     Initially set to \code{NA}.
#' @slot mD Loss rates from pool D (returned to either N or D).
#'     Initially set to \code{NA}.
#' @slot mV Loss rates from pool V (returned to either N or D).
#'     Initially set to \code{0.1}.
#' @slot mH Loss rates from pool H (returned to either N or D).
#'     Initially set to \code{0.1}.
#' @slot mR Loss rates from pool R (returned to either N or D).
#'     Initially set to \code{0.1}.
#' @slot mM Loss rates from pool M (returned to either N or D).
#'     Initially set to \code{0.5}.
#' @slot kP Carrying capacity for P. Initially set to \code{8000}.
#' @slot kV Carrying capacity for V. Only used for model A. Initially set to \code{162}.
#' @slot kH Carrying capacity for H. Only used for model A. Initially set to \code{48}.
#' @slot kR Carrying capacity for R. Only used for model A. Initially set to \code{26}.
#' @slot hP Handing time for P. Only for model B. Initially set to \code{1}.
#' @slot hD Handing time for D. Only for model B. Initially set to \code{1}.
#' @slot hR Handing time for R. Only for model B. Initially set to \code{1}.
#' @slot aNP Uptake rate for NP. Initially set to \code{NA}.
#' @slot aDV Uptake rate for DV. Initially set to \code{NA}.
#' @slot aPH Uptake rate for PH. Initially set to \code{NA}.
#' @slot aR Uptake rate for R. Initially set to \code{NA}.
#' @slot model Which model to use ("A" or "B"). Initially set to \code{"A"}.
#'
#'
#'
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#'
#'
#' @importFrom R6 R6Class
#' @importFrom rootSolve multiroot
#' @importFrom deSolve ode
#' @importFrom dplyr as_tibble
#' @export
#' @keywords data
#'
#'
#' @examples
#' library(tidyverse)
#'
#' # Initialize model (A is default)
#' foodweb_A = web$new(model="A")
#'
#' # Solve for unknown values
#' foodweb_A$eq_solve()
#'
#' # Viewing class:
#' foodweb_A
#'
#' # Outputting class values:
#' foodweb_A$values()
#'
#' # Test midge pulse:
#' plot(foodweb_A$test_midges(1000, a=1000, b=1, r=1, w=400, d=1), type = 'l',
#'      ylab = "iM")
#'
#'
#' # Solve ODEs
#' output_A = foodweb_A$ode_solve(tmax = 1000, a=1000, b=1, r=1, w=400, d=1)  %>%
#'     gather('pool', 'biomass', -time)
#'
#' # Plot absolute biomass
#' output_A  %>%
#'     group_by(pool) %>%
#'     # Define 'minb' to set the minimum value for the y-axis.
#'     # This allows different y-scales for different facets, with the ymin set to 0
#'     mutate(minb = 0) %>%
#'     ggplot(aes(time, biomass)) +
#'     facet_wrap(~pool, scales="free_y") +
#'     # The horizontal lines show the initial states
#'     geom_hline(data = foodweb_A$initial_states,
#'                aes(yintercept=biomass), color="firebrick") +
#'     geom_line(size = 1) +
#'     geom_point(aes(time, minb), shape="") +
#'     theme_classic()
#'
#' # Relative to Equilibrium
#' # Note that this gives weird results when there is no deviation from equilibrium
#' # This is probably due to small numerical errors, but is not a big issue
#' output_A %>%
#'     filter(pool!="M") %>%
#'     group_by(pool) %>%
#'     # Scale the biomass relative to the initial state
#'     mutate(biomass_scale = biomass/biomass[1]) %>%
#'     ggplot(aes(time, biomass_scale)) +
#'     facet_wrap(~pool) +
#'     geom_hline(yintercept = 1, color="firebrick4") +
#'     geom_line(size = 1) +
#'     theme_classic()
#'
#'
#'
#'
#'
web <- R6Class(

    "web",

    # =================================================
    #  Public attributes
    # =================================================

    public = list(

        # Desired equilibrium values
        Neq = 34300, Deq = 308700, Peq = 4300,
        Veq = 81, Heq = 24, Req = 13, Meq = 0,
        # Initial states
        N0 = NULL, D0 = NULL, P0 = NULL, V0 = NULL,
        H0 = NULL, R0 = NULL, M0 = NULL,
        # Inputs
        iN = 1000,
        # Loss rates from systems
        lD = 0.1, lP = 0.1, lV = 0.1, lH = 0.1,
        lR = 0.1, lM = 0.1,
        # Loss rates from pool (returned to either N or D)
        mN = 0.002, mP = NA, mD = NA, mV = 0.1, mH = 0.1,
        mR = 0.1, mM = 0.5,
        # Carying Capacities
        kP = 8000,
        kV = 162, kH = 48, kR = 26,  # only for model A
        # # Handling times
        hP = 1, hD = 1, hR = 1,      # only for model B
        # Uptake rates
        aNP = NA, aDV = NA, aPH = NA, aR = NA,
        # Model A or B
        model = 'A',

        initialize = function(...) {

            # Checking types
            pars <- list(...)
            if (length(pars) > 0) {
                for (i in 1:length(pars)) {
                    if (names(pars)[i] == 'model') {
                        if (! pars[[i]] %in% LETTERS[1:2]) {
                            stop("model must be 'A' or 'B'")
                        }
                    } else if (!is.na(pars[[i]]) & !is.numeric
                               (pars[[i]])) {
                        stop("parameter ", names(pars)[i], " is not",
                             "either NA or numeric, ",
                             "as is required")
                    }
                }
            }
            # Filling pars list with default values
            # (when they're not supplied)
            defaults = private$par_list()
            pars[names(defaults)[!names(defaults) %in% names(pars)]] =
                defaults[!names(defaults) %in% names(pars)]
            # If an initial state is NULL,
            # use that pool's equilibrium value
            sv = private$pool_names()
            for (s in sv) {
                if (is.null(pars[[paste0(s, '0')]])) {
                    pars[[paste0(s, '0')]] = pars[[paste0(s, 'eq')]]
                }
            }

            private$assign_list(pars)
        },

        # Use estimated equilibrium biomass
        # and selected values for certain
        # parameters to solve for unknown values
        # "Known" parameters must be selected with care,
        # to ensure that an equilibrium solution can
        # actually be reached
        eq_solve = function(solve_pars = c("mP", "mD", "aNP", "aDV", "aPH", "aR"),
                            initial_vals = rep(0.1, 6)) {

            stopifnot(length(solve_pars) == 6,
                      length(initial_vals) == 6)
            if (! all(solve_pars %in% private$par_names())) {
                stop("One or more solve parameter names is not",
                     "present in this class.")
            }

            cat("Solving for the following parameters:",
                paste(solve_pars, collapse = ", "), "\n")

            pars = private$par_list()
            pars[solve_pars] = NA

            # Checking for too many unknown values
            if (length(pars[is.na(pars)]) != 6) {
                stop("cannot solve for unknown parameters when the ",
                     "number of unknown ",
                     "(i.e., NA) parameters is != 6; ",
                     "the following are NA: ",
                     paste(names(pars)[is.na(pars)], collapse = ", "))
            }

            start = initial_vals
            names(start) = solve_pars

            par_solve = multiroot(f = private$equil_solve,
                                  start = start, parms = pars,
                                  maxiter = 100,
                                  rtol = 1e-6, atol = 1e-8,
                                  ctol = 1e-8,
                                  useFortran = TRUE,
                                  positive = TRUE,
                                  jacfunc = NULL,
                                  jactype = "fullint",
                                  verbose = FALSE, bandup = 1,
                                  banddown = 1)

            pars[names(par_solve$root)] = as.list(par_solve$root)

            private$assign_list(pars)

            invisible(NULL)
        },


        # Function to Solve ODE
        # tmin and tmax specify the duration over which to
        # run the model
        # tstep specifies the step size
        # a, b, r, w, and d define the midge pulse (see `vignette("smooth_pulse", "mtf")`)
        ode_solve = function(tmax, a, b, r, w, d, tstep = 1){

            # pars and init give the parameter values and
            # initial states
            pars <- private$par_list()
            pars$midges <- function(t_) private$midge_pulse(t_, a, b, r, w, d)

            # Checking for NA values
            if (length(pars[is.na(pars)]) > 0) {
                stop(sprintf("You still have the following NA parameters: %s.\n%s",
                     paste(names(pars)[is.na(pars)], collapse = ", "),
                     "Please solve for these using web$eq_solve() first."))
            }

            init = unlist(pars[paste0(private$pool_names(), '0')])
            names(init) = private$pool_names()

            solved_ode = ode(init, seq(0, tmax, tstep),
                                      private$diff_eq, pars)
            solved_ode = as_tibble(as.data.frame(solved_ode))

            return(solved_ode)

        },


        values = function() {
            pn = private$par_list()
            return(pn)
        },


        test_midges = function(tmax, a, b, r, w, d, tstep = 1) {
            return(private$midge_pulse(seq(1, tmax, tstep), a, b, r, w, d))
        },


        print = function(...) {
            cat('<< Class web >>\n')
            cat('< Foodweb model', self$model, '>\n')
            cat('Nutrient pools: \n')
            cat('  - N: soil (plant-available) \n')
            cat('  - D: detritus \n')
            cat('  - P: plants \n')
            cat('  - V: detritivores \n')
            cat('  - H: herbivores \n')
            cat('  - R: predators \n')
            cat('  - M: midges \n')

            cat('\nParameters used in model:\n')
            cat('  - X0: starting value for pool X\n')
            cat('  - Xeq: desired equilibrium value for pool X\n')
            cat('  - iN: input to soil nutrient pool\n')
            cat('  - lX: loss rates from system from pool X
                (X: D,P,V,H,R,M)\n')
            cat('  - mX: loss rates from pool X, returned to either
                N or D\n',
                '   (X: N,D,P,V,H,R,M)\n')
            if (self$model == 'A') {
                cat('  - kX: carrying capacities for pool X
                    (X: P,V,H,R)\n')
            } else {
                cat('  - kX: carrying capacities for pool X (X: P)\n')
                cat('  - hX: handling time for pool X (X: P,D,R)\n')

            }
            cat('  - aNP: maximum uptake rate from pool X to pool Y
                (XY: NP, CV, PH)\n')
            cat('  - aR: maximum uptake rate to predators;
                same value for V,H,M\n')
        }

    ),

    # =================================================
    #  Private attributes
    # =================================================

    private = list(

        # Define function for steady states to solve
        equil_solve = function(start, all_params){

            all_params = as.list(all_params)

            # Set starting values
            all_params[names(start)] = as.list(start)

            # Pass parameters to dynamic equations
            if (self$model == 'A') {
                output = with(all_params,
                              c(sN = iN - aNP*Neq*Peq*(1-Peq/kP) + (1-lD)*mD*Deq - mN*Neq,
                                sD = (1-lP)*mP*Peq + (1-lV)*mV*Veq + (1-lH)*mH*Heq +
                                    (1-lR)*mR*Req - aDV*Deq*Veq*(1-Veq/kV) - mD*Deq,
                                sP = aNP*Neq*Peq*(1-Peq/kP) - aPH*Peq*Heq*(1-Heq/kH) - mP*Peq,
                                sV = aDV*Deq*Veq*(1-Veq/kV) - (aR*Veq*Req)*(1-Req/kR) - mV*Veq,
                                sH = aPH*Peq*Heq*(1-Heq/kH) - (aR*Heq*Req)*(1-Req/kR) - mH*Heq,
                                sR = (aR*Veq*Req + aR*Heq*Req)*(1-Req/kR) - mR*Req
                              ))
            } else {
                output = with(all_params,
                              c(sN = iN - aNP*Neq*Peq*(1-Peq/kP) + (1-lD)*mD*Deq - mN*Neq,
                                sD = (1-lP)*mP*Peq + (1-lV)*mV*Veq + (1-lH)*mH*Heq +
                                    (1-lR)*mR*Req - aDV*Deq*Veq/(1+aDV*hD*Deq) - mD*Deq,
                                sP = aNP*Neq*Peq*(1-Peq/kP) - aPH*Peq*Heq/(1+aPH*hP*Peq) - mP*Peq,
                                sV = aDV*Deq*Veq/(1+aDV*hD*Deq) - (aR*Veq*Req)/(1+aR*hR*(Veq+Heq)) -
                                    mV*Veq,
                                sH = aPH*Peq*Heq/(1+aPH*hP*Peq) - (aR*Heq*Req)/(1+aR*hR*(Veq+Heq)) -
                                    mH*Heq,
                                sR = (aR*Veq*Req + aR*Heq*Req)/(1+aR*hR*(Veq+Heq)) - mR*Req
                              ))
            }

            return(output)

        },

        # State variable parameter names
        pool_names = function() {
            c(
                "N",
                "D",
                "P",
                "V",
                "H",
                "R",
                "M"
            )
        },

        # All variable names
        par_names = function() {
            c(
                "model",
                "Neq",
                "Deq",
                "Peq",
                "Veq",
                "Heq",
                "Req",
                "Meq",
                "N0",
                "D0",
                "P0",
                "V0",
                "H0",
                "R0",
                "M0",
                "iN",
                "lD",
                "lP",
                "lV",
                "lH",
                "lR",
                "lM",
                "mN",
                "mP",
                "mD",
                "mV",
                "mH",
                "mR",
                "mM",
                "kP",
                "kV",    # only for model A
                "kH",    # only for model A
                "kR",    # only for model A
                "hP",    # only for model B
                "hD",    # only for model B
                "hR",    # only for model B
                "aNP",
                "aDV",
                "aPH",
                "aR"
            )
        },

        # Return all parameters as a single list
        par_list = function() {
            L = list()
            pn = private$par_names()
            for (p in pn) {
                L[[p]] = self[[p]]
            }
            return(L)
        },

        # Assign all parameters from a single list
        assign_list = function(L) {

            pn = private$par_names()
            for (p in pn) {
                self[[p]] = L[[p]]
            }
        },

        # Define differential equations
        diff_eq = function(t, y, pars) {

            iM = pars$midges(t)

            if (self$model == 'A') {
                output =
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
            } else {
                output =
                    with(as.list(pars),
                         with(as.list(y),
                              {
                                  c(N = iN - aNP*N*P*(1-P/kP) + (1-lD)*mD*D - mN*N,
                                    D = (1-lP)*mP*P + (1-lV)*mV*V + (1-lH)*mH*H +
                                        (1-lR)*mR*R + (1-lM)*mM*M -
                                        aDV*D*V/(1+aDV*hD*D) - mD*D,
                                    P = aNP*N*P*(1-P/kP) - aPH*P*H/(1+aPH*hP*P) - mP*P,
                                    V = aDV*D*V/(1+aDV*hD*D) -
                                        (aR*V*R)/(1+aR*hR*(V+H+M)) - mV*V,
                                    H = aPH*P*H/(1+aPH*hP*P) -
                                        (aR*H*R)/(1+aR*hR*(V+H+M)) - mH*H,
                                    R = (aR*V*R + aR*H*R + aR*M*R)/(1+aR*hR*(V+H+M)) -
                                        mR*R,
                                    M = iM - mM*M - (aR*M*R)/(1+aR*hR*(V+H+M)))
                              }))
            }

            return(list(output))
        },

        midge_pulse = function(t, a, b, r, w, d) {
            max_adj <- sin(2 * pi * (0.75 - 1 / (3 * r)))
            max <- b * (1 + exp(-a * (1 + max_adj)))
            sines <- sin(2 * pi * (8 * t - 8 * d * w + 9 * r * w)/(12 * r * w))
            f <- max / (1+exp(a * (sines - max_adj)))
            return(f)
        }

    ),

    active = list(

        initial_states = function(){
            data.frame(pool = c("N","D","V","P","H","R","M"),
                       biomass =  c(self$N0,self$D0,self$V0,self$P0,self$H0,self$R0,self$M0))

        }

    )
)




multi_web <- function(par_list, expand = FALSE) {

    if (!expand) {
        if (diff(range(sapply(par_list, length))) != 0) {
            stop("input parameter list must have vectors of the same length",
                 "if `expand == FALSE`")
        }
        web_list <- lapply(1:length(par_list[[1]]),
                           function(i) {
                               input_list <- lapply(par_list, function(x) x[[i]])
                               new_web <- do.call(web$new, input_list)
                               return(new_web)
                           })
    } else {

        par_df <- do.call(expand.grid, par_list)


        web_list <- lapply(
            1:nrow(par_df),
            function(i) {
                input_list <- as.list(par_df[i,])
                attr(input_list, 'out.attrs') <- NULL
                new_web <- do.call(web$new, input_list)
                return(new_web)
            })
    }

    return(web_list)
}

