library(rootSolve)


diff_eq <- function(x, parms) {

    parms[names(x)] <- as.list(x)

    # output <- with(as.list(parms),
    #                     c(sN = - mN*N_t - ((aNP*N_t*P_t)/(1+bP*P_t)) + aND*D_t*V_t + G,
    #                       sD = mV*V_t + mP*P_t + mH*H_t + mR*R_t - aND*D_t*V_t - aVD*V_t*D_t + f*M,
    #                       sV = - mV*V_t -((aVR*V_t*R_t)/(1 + bR*R_t)) + aVD*V_t*D_t,
    #                       sP = - mP*P_t + ((aNP*N_t*P_t)/(1 + bP*P_t)) - aPH*P_t*H_t,
    #                       sH = - mH*H_t - ((aHR*H_t*R_t)/(1 + bR*R_t)) + aPH*P_t*H_t,
    #                       sR = - mR*R_t + ((aHR*H_t*R_t)/(1 + bR*R_t)) +
    #                           ((aVR*V_t*R_t)/(1 + bR*R_t)) + (1 - f)*M)
    # )
    output <- with(as.list(parms),
                        c(sN = - mN*N_t - ((aNP*N_t*P_t)/(1+bP*P_t)) + aND*D_t*V_t + G,
                          sD = mV*V_t + mP*P_t + mH*H_t + mR*R_t - aND*D_t*V_t - aVD*V_t*D_t, # + f*M,
                          sV = - mV -((aVR*R_t)/(1 + bR*R_t)) + aVD*D_t,
                          sP = - mP + ((aNP*N_t)/(1 + bP*P_t)) - aPH*H_t,
                          sH = - mH - ((aHR*R_t)/(1 + bR*R_t)) + aPH*P_t,
                          sR = - mR + ((aHR*H_t)/(1 + bR*R_t)) +
                              ((aVR*V_t)/(1 + bR*R_t)) )#+ (1 - f)*M)
    )

    return(output)
}

# Initial N pools for variables
N_t <- 343000 #nutrient pool
D_t <- 114000 #detritus
V_t <- 81 #detritivores
P_t <- 4300 #plants
H_t <- 24 #herbivores
R_t <- 13 #predators

parms <- list(mN = 100/N0, mV = 2/V0, mP = 0.1/P0, mH = 2/H0, mR = 2/R0,
              G = 100, bP = 100, bR = 0.01, M = 0, f = 0 * 2/3,
              N_t = 343000, D_t = 114000, V_t = 81, P_t = 4300, H_t = 24, R_t = 13,
              aND = NA, aVD = NA, aNP = NA, aPH = NA, aVR = NA, aHR = NA)

start <- c(1.082954e-07, 2.165909e-07, 2.915452e-06, 2.325581e-04, 1.234568e-02, 4.166667e-02)
names(start) <- c("aND", "aVD", "aNP", "aPH", "aVR", "aHR")

pmatrix; unlist(parms)

# start <- c(aND = 2.27e-7,
#            aVD = 4.33e-7,
#            aNP = 6.12e-4,
#            aPH = 1.94e-5,
#            aVR = 2.15e-3,
#            aHR = 1.13e-7)
diff_eq(start, parms)


z <- multiroot(f = diff_eq, start = start, parms = parms, maxiter = 100, rtol = 1e-6,
               atol = 1e-8, ctol = 1e-8, useFortran = TRUE, positive = T, jacfunc = NULL,
               jactype = "fullint", verbose = FALSE, bandup = 1, banddown = 1)
z

