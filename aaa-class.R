
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
        # Midge function 
        #(have to wrap inside a list bc it's unchangeable otherwise)
        iM_func = list(function(t) 0),
        # Model A or B
        model = 'A',
        
        initialize = function(..., do_solve = TRUE, 
                              initial_vals = rep(0.1, 6)) {
            
            # Checking types
            pars <- list(...)
            if (length(pars) > 0) {
                for (i in 1:length(pars)) {
                    if (names(pars)[i] == 'iM_func') {
                        if (!is.function(pars[[i]])) {
                            stop("parameter iM_func must be a
                                 function")
                        } else if (!identical(formalArgs(pars[[i]]),
                                              't')) {
                            stop("parameter iM_func must be a",
                                 "function with the only ",
                                 "argument being 't'")
                        }
                    } else if (names(pars)[i] == 'model') {
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
            
            # Checking for too many unknown values
            if (length(pars[is.na(pars)]) != 6 & do_solve) {
                stop("cannot solve for unknown parameters when the",
                     "number of unknown ",
                     "(i.e., NA) parameters is != 6; ",
                     "by default the following are NA: ",
                     "mP, mD, aNP, aDV, aPH, and aR")
            }
            if (length(pars[is.na(pars)]) > 0 & !do_solve) {
                stop("if not solving for unknown parameters, you",
                     "cannot have any ",
                     "parameters as NA; ",
                     "by default the following are NA: ",
                     "mP, mD, aNP, aDV, aPH, and aR")
            }
            
            # Use estimated equilibrium biomass 
            # and selected values for certain 
            # parameters to solve for unknown values
            # "Known" parameters must be selected with care, 
            # to ensure that an equilibrium solution can 
            # actually be reached
            if (do_solve) {
                
                start = initial_vals
                names(start) = names(pars)[is.na(pars)]
                
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
            }
            
            private$assign_list(pars)
        },
        
        
        re_solve = function(solve_pars = c("mP", "mD", "aNP", "aDV",
                                           "aPH", "aR"),
                            initial_vals = rep(0.1, 6)) {
            
            stopifnot(length(solve_pars) == 6, 
                      length(initial_vals) == 6)
            if (! all(solve_pars %in% private$par_names())) {
                stop("One or more solve parameter names is not",
                     "present in this class.")
            }
            
            pars = private$par_list()
            pars[solve_pars] = NA
            
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
        ode_solve = function(tmin, tmax, tstep){
            
            # parms and init give the parameter values and 
            # initial states
            parms = private$par_list()
            
            init = unlist(parms[paste0(private$pool_names(), '0')])
            names(init) = private$pool_names()
            
            solved_ode = ode(init, seq(tmin,tmax,tstep), 
                             private$diff_eq, parms)
            solved_ode = as_tibble(as.data.frame(solved_ode))
            
            return(solved_ode)
            
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
            cat('  - iM_func: function for midge input at time t\n')
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
                "aR",
                "iM_func",
                "model"
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
            
            iM = parms$iM_func(t)
            
            if (self$model == 'A') {
                output = 
                    with(as.list(parms), 
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
                    with(as.list(parms),
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
        
        par_df <- do.call(expand.grid, par_list[names(par_list) != 'iM_func'])
        if ('iM_func' %in% names(par_list)) {
            par_df_list <- lapply(1:length(par_list[['iM_func']]), 
                   function(i) {
                       par_df %>% 
                           as_tibble %>% 
                           mutate(iM_func = par_list[['iM_func']][i])
                   })
            par_df <- bind_rows(par_df_list)
        }
        
        
        web_list <- lapply(
            1:nrow(par_df), 
            function(i) {
                input_list <- as.list(par_df[i,])
                attr(input_list, 'out.attrs') <- NULL
                if ('iM_func' %in% names(input_list)) {
                    input_list[['iM_func']] <- input_list[['iM_func']][[1]]
                }
                new_web <- do.call(web$new, input_list)
                return(new_web)
            })
    }
    
    return(web_list)
}

