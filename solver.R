library(rootSolve)



diff_eq <- function(y, parms) {
    
    output <- with(as.list(parms),
                   with(as.list(y),
                        c(N_t = - mN*N_t - ((aNP*N_t*P_t)/(1+bP*P_t)) + aND*D_t*V_t + G,
                          D_t = mV*V_t + mP*P_t + mH*H_t + mR*R_t - aND*D_t*V_t - 
                              aVD*V_t*D_t + f*M,
                          V_t = - mV*V_t -((aVR*D_t*R_t)/(1 + bR*R_t)) + aVD*V_t*D_t,
                          P_t = - mP*P_t + ((aNP*N_t*P_t)/(1 + bP*P_t)) - aPH*P_t*H_t,
                          H_t = - mH*H_t - ((aHR*H_t*R_t)/(1 + bR*R_t)) + aPH*P_t*H_t,
                          R_t = - mR*R_t + ((aHR*H_t*R_t)/(1 + bR*R_t)) + 
                              ((aVR*D_t*R_t)/(1 + bR*R_t)) + (1 - f)*M)
                   )
    )
    
    return(output)
}

# Initial N pools for variables
N0 <- 343000 #nutrient pool
D0 <- 114000 #detritus
V0<- 81 #detritivores
P0<- 4300 #plants
H0 <- 24 #herbivores
R0 <- 13 #predators

parms <- list(mN = 100/N0, mV = 2/V0, mP = 0.1/P0, mH = 0.001 * 2/H0, mR = 0.1 * 2/R0,
              G = 10, bP = 100, bR = 0, M = 0, f = 2/3,
              aND = NA * 2.27e-7, 
              aVD = NA * 4.33e-7, 
              aNP = NA * 1e-3 * 6.12e-4,
              aPH = NA * 1e-2 * 1.94e-5, 
              aVR = NA * 1e-4 * 2.15e-3, 
              aHR = NA * 1.13e-7)

start <- c(N_t = 343000, D_t = 114000, V_t = 81, P_t = 4300, H_t = 24, R_t = 13)


z <- multiroot(f = diff_eq, start = start, parms = parms, maxiter = 100, rtol = 1e-6, 
               atol = 1e-8, ctol = 1e-8, useFortran = TRUE, positive = T,jacfunc = NULL,
               jactype = "fullint", verbose = FALSE, bandup = 1, banddown = 1)
