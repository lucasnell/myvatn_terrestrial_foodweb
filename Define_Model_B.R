
# ======================================
# Solve for unknown parameters
# ======================================
# Use estimated equilibrium biomass and selected values for certain parameters to solve 
# for unknown values
# "Known" parameters must be selected with care, to ensure that an equilibrium solution 
# can actually be reached

# Define function for steady states to solve
equil_solve = function(start, all_params){
  
  # Set starting values
  all_params[names(start)] = as.list(start)
  
  # Pass parameters to dynamic equations
  output = as.list(all_params) %>% with(
    c(sN = iN - aNP*N*P*(1-P/kP) + (1-lD)*mD*D,
      sD = (1-lP)*mP*P + (1-lV)*mV*V + (1-lH)*mH*H + (1-lR)*mR*R - 
        aDV*D*V/(1+aDV*hD*D) - mD*D, 
      sP = aNP*N*P*(1-P/kP) - aPH*P*H/(1+aPH*hP*P) - mP*P,
      sV = aDV*D*V/(1+aDV*hD*D) - (aR*V*R)/(1+aR*hR*(V+H)) - mV*V,
      sH = aPH*P*H/(1+aPH*hP*P) - (aR*H*R)/(1+aR*hR*(V+H)) - mH*H,
      sR = (aR*V*R + aR*H*R)/(1+aR*hR*(V+H)) - mR*R
    ))
  
  return(output)
  
}


# Define equilibrium states
# Define known and unknown parameters
# Merge
equil_states = list(N = 343000, D = 114000, P = 4300, V = 81, H = 24,  R = 13)

params <- list(
  # Inputs 
  iN = 10, 
  # Loss rates from systems
  lD = 0.1, lP = 0.1, lV = 0.1, lH = 0.1, lR = 0.1,
  # Loss rates from pool (returned to either N or D)
  mP=NA, mD=NA, mV=0.1, mH=0.1, mR=0.1,
  # Carying Capacities
  kP=8000, 
  # Handling Times
  hP = 1, hD = 1, hR = 1,
  # Uptake rates
  aNP=NA, aDV=NA, aPH=NA, aR=NA)

all_params = c(equil_states,params)

# Initial values for unknown rates
start = rep(0.01,6)
names(start) = c("mD","mP","aNP","aDV","aPH","aR")

# Test to make sure function returns values when provided all_params
equil_solve(start, all_params)

# Solve for unknown parameters
z = multiroot(f = equil_solve, start = start, parms = all_params, maxiter = 100, 
              rtol = 1e-6, atol = 1e-8, ctol = 1e-8, useFortran = TRUE, 
              positive = T, jacfunc = NULL,jactype = "fullint", 
              verbose = FALSE, bandup = 1, banddown = 1)

# Examine solution
round(z$root,10)

# Add solved values to final parameter set
params_solve = params
params_solve[names(z$root)] = as.list(z$root)





# ======================================
# Define Functions to Run Model
# ======================================

# Create function for time-varying midge pulse
iM_func = function(t,pulse,pulse_tmin,pulse_tmax) 
{ifelse(t>pulse_tmin&t<pulse_tmax,pulse,0)}
params_solve$iM_func = iM_func

# Define differential equations
diff_eq = function(t, y, parms) {
  
  iM = parms %>% with(parms[["iM_func"]](t,pulse,pulse_tmin,pulse_tmax))
  
  output = as.list(parms) %>% 
    with(as.list(y) %>% 
           with({
             c(N = iN - aNP*N*P*(1-P/kP) + (1-lD)*mD*D,
               D = (1-lP)*mP*P + (1-lV)*mV*V + (1-lH)*mH*H + (1-lR)*mR*R + 
                 (1-lM)*mM*M - aDV*D*V/(1+aDV*hD*D) - mD*D, 
               P = aNP*N*P*(1-P/kP) - aPH*P*H/(1+aPH*hP*P) - mP*P,
               V = aDV*D*V/(1+aDV*hD*D) - (aR*V*R)/(1+aR*hR*(V+H+M)) - mV*V,
               H = aPH*P*H/(1+aPH*hP*P) - (aR*H*R)/(1+aR*hR*(V+H+M)) - mH*H,
               R = (aR*V*R + aR*H*R + aR*M*R)/(1+aR*hR*(V+H+M)) - mR*R,
               M = iM - mM*M - (aR*M*R)/(1+aR*hR*(V+H+M))
             )}))
  
  return(list(output))
}

# Function to Solve ODE
ode_solve = 
  function(tmin,tmax,tstep,pulse,pulse_tmin,pulse_tmax,parms,init){
    parms$pulse = pulse
    parms$pulse_tmin = pulse_tmin
    parms$pulse_tmax = pulse_tmax
    ss = ode(init, seq(tmin,tmax,tstep), diff_eq, parms)
    
    return(
      if(tmax>=pulse_tmax&tmin<=pulse_tmin) ss else 
        warning("Time Mismatch"))
  }








