
library(R6)

web <- R6Class(
    
    "web",
    
    # =================================================
    #  Public attributes
    # =================================================
    
    public = list(
        
        # Initial states
        N = 343000, D = 114000, P = 4300, V = 81, H = 24,  R = 13, M = 0,
        # Inputs
        iN = 1000,
        # Loss rates from systems
        lD = 0.1, lP = 0.1, lV = 0.1, lH = 0.1, lR = 0.1, lM = 0.1,
        # Loss rates from pool (returned to either N or D)
        mP = NA, mD = NA, mV = 0.1, mH = 0.1, mR = 0.1, mM = 0.5,
        # Carying Capacities
        kP = 8000, 
        kV = 162, kH = 48, kR = 26,  # comment for model B
        # # Handling times
        # hP = 1, hD = 1, hR = 1,    # uncomment for model B
        # Uptake rates
        aNP = NA, aDV = NA, aPH = NA, aR = NA,
        # Midge function (have to wrap inside a list bc it's unchangeable otherwise)
        iM_func = list(function(t) 0),
        
        initialize = function(..., do_solve = TRUE, initial_vals = rep(0.1, 6)) {
            
            # Checking types
            pars <- list(...)
            if (length(pars) > 0) {
                for (i in 1:length(pars)) {
                    if (names(pars)[i] == 'iM_func') {
                        if (!is.function(pars[[i]])) {
                            stop("parameter iM_func must be a function")
                        } else if (!identical(formalArgs(pars[[i]]), 't')) {
                            stop("parameter iM_func must be a function with the only ",
                                 "argument being 't'")
                        }
                    } else if (!is.na(pars[[i]]) & !is.numeric(pars[[i]])) {
                        stop("parameter ", names(pars)[i], " is not either NA or numeric, ",
                             "as is required")
                    }
                }
            }
            # Filling pars list with default values when they're not supplied
            defaults = private$par_list()
            pars[names(defaults)[!names(defaults) %in% names(pars)]] = 
                defaults[!names(defaults) %in% names(pars)]
            
            # Checking for too many unknown values
            if (length(pars[is.na(pars)]) != 6 & do_solve) {
                stop("cannot solve for unknown parameters when the number of unknown ",
                     "(i.e., NA) parameters is != 6; ",
                     "by default the following are NA: ",
                     "mP, mD, aNP, aDV, aPH, and aR")
            }
            if (length(pars[is.na(pars)]) > 0 & !do_solve) {
                stop("if not solving for unknown parameters, you cannot have any ",
                     "parameters as NA; ",
                     "by default the following are NA: ",
                     "mP, mD, aNP, aDV, aPH, and aR")
            }
            
            if (do_solve) {
                start = initial_vals
                names(start) = names(pars)[is.na(pars)]
                
                par_solve = multiroot(f = private$equil_solve, 
                                      start = start, parms = pars,
                                      maxiter = 100,
                                      rtol = 1e-6, atol = 1e-8, ctol = 1e-8,
                                      useFortran = TRUE, positive = TRUE, jacfunc = NULL,
                                      jactype = "fullint", verbose = FALSE, bandup = 1,
                                      banddown = 1)
                
                pars[names(par_solve$root)] = as.list(par_solve$root)
            }
            
            private$assign_list(pars)
        },
        
        
        re_solve = function(solve_pars = c("mP", "mD", "aNP", "aDV", "aPH", "aR"),
                            initial_vals = rep(0.1, 6)) {
            
            stopifnot(length(solve_pars) == 6, length(initial_vals) == 6)
            if (! all(solve_pars %in% private$par_names())) {
                stop("One or more solve parameter names is not present in this class.")
            }
            
            pars = private$par_list()
            pars[solve_pars] = NA
            
            start = initial_vals
            names(start) = solve_pars
            
            par_solve = multiroot(f = private$equil_solve, 
                                  start = start, parms = pars,
                                  maxiter = 100,
                                  rtol = 1e-6, atol = 1e-8, ctol = 1e-8,
                                  useFortran = TRUE, positive = TRUE, jacfunc = NULL,
                                  jactype = "fullint", verbose = FALSE, bandup = 1,
                                  banddown = 1)
            
            pars[names(par_solve$root)] = as.list(par_solve$root)
            
            private$assign_list(pars)
            
            invisible(NULL)
        },
        
        
        # Function to Solve ODE
        ode_solve = function(tmin, tmax, tstep){
            
            parms = private$par_list()
            
            init = unlist(parms[c("N", "D", "P", "V", "H", "R", "M")])
            
            solved_ode = ode(init, seq(tmin,tmax,tstep), private$diff_eq, parms)
            solved_ode = as_tibble(as.data.frame(solved_ode))
            
            return(solved_ode)
            
        }  #,
        
        
        # print = function(...) {
        #     cat('Yup, foodweb from Myvatn.\n')
        # }
        
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
            output = with(all_params,
                          c(sN = iN - aNP*N*P*(1-P/kP) + (1-lD)*mD*D,
                            sD = (1-lP)*mP*P + (1-lV)*mV*V + (1-lH)*mH*H + (1-lR)*mR*R - 
                                aDV*D*V*(1-V/kV) - mD*D, 
                            sP = aNP*N*P*(1-P/kP) - aPH*P*H*(1-H/kH) - mP*P,
                            sV = aDV*D*V*(1-V/kV) - (aR*V*R)*(1-R/kR) - mV*V,
                            sH = aPH*P*H*(1-H/kH) - (aR*H*R)*(1-R/kR) - mH*H,
                            sR = (aR*V*R + aR*H*R)*(1-R/kR) - mR*R
                          ))
            
            return(output)
            
        },
        
        
        par_names = function() {
            c(
                "N",
                "D",
                "P",
                "V",
                "H",
                "R",
                "M",
                "iN",
                "lD",
                "lP",
                "lV",
                "lH",
                "lR",
                "lM",
                "mP",
                "mD",
                "mV",
                "mH",
                "mR",
                "mM",
                "kP",
                "kV",    # comment for model B
                "kH",    # comment for model B
                "kR",    # comment for model B
                # "hP",  # uncomment for model B
                # "hD",  # uncomment for model B
                # "hR",  # uncomment for model B
                "aNP",
                "aDV",
                "aPH",
                "aR",
                "iM_func"
            )
        },
        
        # Return all parameters as a single list
        par_list = function() {
            L = list()
            pn = private$par_names()
            for (p in pn) {
                if (p != 'iM_func') {
                    L[[p]] = self[[p]]
                } else {
                    L[[p]] = self[[p]][[1]]
                }
            }
            return(L)
        },
        
        # Assign all parameters from a single list
        assign_list = function(L) {
            
            pn = private$par_names()
            for (p in pn) {
                if (p != 'iM_func') {
                    self[[p]] = L[[p]]
                } else {
                    self[[p]] = list(L[[p]])
                }
            }
        },
        
        # Define differential equations
        diff_eq = function(t, y, parms) {
            
            iM = with(parms, iM_func(t))
            
            output = 
                with(as.list(parms), 
                     with(as.list(y), 
                          {
                              c(N = iN - aNP*N*P*(1-P/kP) + (1-lD)*mD*D,
                                D = (1-lP)*mP*P + (1-lV)*mV*V + (1-lH)*mH*H + (1-lR)*mR*R + 
                                    (1 - lM)*mM*M - aDV*D*V*(1-V/kV) - mD*D, 
                                P = aNP*N*P*(1-P/kP) - aPH*P*H*(1-H/kH) - mP*P,
                                V = aDV*D*V*(1-V/kV) - (aR*V*R)*(1-R/kR) - mV*V,
                                H = aPH*P*H*(1-H/kH) - (aR*H*R)*(1-R/kR) - mH*H,
                                R = (aR*V*R + aR*H*R + 0.5*aR*M*R)*(1-R/kR) - mR*R,
                                M = iM - mM*M - (0.5*aR*M*R)*(1-R/kR))
                          }))
            
            return(list(output))
        }
        
    )
)
