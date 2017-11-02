#*************************************************
##McCary et al. Nov 2017
##Food-web model for Myvatn
#*************************************************

library(ggplot2)
library(dplyr)
library(tidyr)
library(deSolve)

##The model tracks the change of 6 components of a food web through time.
##At time t, the total nitrogen pool is denoted by N(t),
##the total amount of detritus is denoted by D(t),
##the abundance of detritivores by V(t), plants by P(t), herbivores by H(t),
##and the predators by P(t).

##We assume that at each time step, except for the detritus,
##each component has specific mortality rate denoted by mX where X indicate
##with component it is. For instance, mD denotes the per unit mortality rate
##of the detritivores. A consumed unit of a component X is converted into a
##consumer component Y at a rate aXY. For instance, one unit of nutrient is
##converted into aND unit of detritus.

##Finally, we assume that consumptions by plants and predators follow a
##Holling type II functional response with an attack rate of bP and bR, respectively.
##At each time step, we add nitrogen and midges denoted by G and M, respectively.
##A fraction f of the midge goes to the detritivores whereas the remaining fraction (1-f)
##goes to the plant. Under these assumptions, the "abundances" at time t +1 are

##Define parameters
#Initial N pools for variables
N0 <- 343000 #nutrient pool
D0 <- 114000 #detritus
V0<- 81 #detritivores
P0<- 4300 #plants
H0 <- 24 #herbivores
R0 <- 13 #predators

##mortality or rates of loss [Detritus doesnt die!]
mN <- 100/N0 #nutrient loss
mV <- 2/V0 #detritivore mortality
mP <- 0.1/P0 #plant mortality
mH <- 2/H0 #herbivore mortality
mR <- 2/R0 #predator mortality


## Unknown transitions rates
aND <- rnorm(1)  # NA
aVD <- rnorm(1)  # NA
aNP <- rnorm(1)  # NA
aPH <- rnorm(1)  # NA
aVR <- rnorm(1)  # NA
aHR <- rnorm(1)  # NA


##Density dependence
bP <- 100 #Plants
bR <- 0 #Predators
G <- 100 #Nutrients added
M <- 0 #Midge deposition

##Midge fraction
f <- 2/3

##Max time
nt <- 100


sim <- function(aND, aVD, aNP, aPH, aVR, aHR, N0 = 343000, D0 = 114000, V0 = 81,
    P0 = 4300, H0 = 24, R0 = 13, mN = NULL, mV = NULL, mP = NULL, mH = NULL, mR = NULL,
    bP = 100, bR = 0, G = 100, M = 0, f = 2/3, nt = 100) {
    if(is.null(mN)) mN <- 100/N0
    if(is.null(mV)) mV <- 2/V0
    if(is.null(mP)) mP <- 0.1/P0
    if(is.null(mH)) mH <- 2/H0
    if(is.null(mR)) mR <- 2/R0

    #Create data frame to store results
    data = data.frame(
    N = c(N0, rep(0, nt-1)),  # nutrients
    D = D0,  # detritus
    V = V0,  # detritivores
    P = P0,  # plants
    H = H0,  # herbivores
    R = R0)  # predators

    for (t in 1:(nt-1)){
        data$N[t+1] <- data$N[t] - mN*data$N[t] - ((aNP*data$N[t]*data$P[t])/(1+bP*data$P[t])) +
            aND*data$D[t]*data$V[t] + G
        data$D[t+1] <- data$D[t] + mV*data$V[t] + mP*data$P[t] + mH*data$H[t] +
            mR*data$R[t] - aND*data$D[t]*data$V[t] - aVD*data$V[t]*data$D[t] + f*M
        data$V[t+1] <- data$V[t] - mV*data$V[t] -((aVR*data$D[t]*data$R[t])/(1 +
            bR*data$R[t])) + aVD*data$V[t]*data$D[t]
        data$P[t+1] <- data$P[t] - mP*data$P[t] + ((aNP*data$N[t]*data$P[t])/(1 +
            bP*data$P[t])) + aPH*data$P[t]*data$H[t]
        data$H[t+1] <- data$H[t] + mH*data$H[t] - ((aHR*data$H[t]*data$R[t])/(1 +
            bR*data$R[t])) + aPH*data$P[t]*data$H[t]
        data$R[t+1] <- data$R[t] - mR*data$R[t] + ((aHR*data$H[t]*data$R[t])/(1 +
            bR*data$R[t])) + ((aVR*data$D[t]*data$R[t])/(1 + bR*data$R[t])) + (1 - f)*M
    }

    return(data)
}


aND <- 2.27e-7
aVD <- 4.33e-7
aNP <- 6.12e-4
aPH <- 1.94e-5
aVR <- 2.15e-3
aHR <- 1.13e-7

matplot(1:nt, sim(aND, aVD, aNP, aPH, aVR, aHR), type = 'l')




nt <- 100

sim(aND, aVD, aNP, aPH, aVR, aHR, M = 0, G = 0, nt = 15)

lapply(c(0, N0 / 1e6, N0 / 1e5), function(x) sim(aND, aVD, aNP, aPH, aVR, aHR, M = x, nt = nt)) %>%
    lapply(function(x) mutate(x, time = 1:nt)) %>%
    bind_rows %>%
    mutate(M = factor(rep(0:2, each = nt), levels = 0:2)) %>%
    gather('pool', 'biomass', -time, -M) %>%
    as_tibble %>%
    ggplot(aes(time, log(biomass), color = pool)) +
    geom_line(size = 1) +
    theme_classic() +
    facet_wrap(~ M) +
    scale_color_brewer(palette = 'Accent')


diff_eq <- function(t, y, parms) {

    M = parms[["M_func"]](t)
    output <- with(as.list(parms),
        with(as.list(y),
            {
            c(N_t = - mN*N_t - ((aNP*N_t*P_t)/(1+bP*P_t)) + aND*D_t*V_t + G,
            D_t = mV*V_t + mP*P_t + mH*H_t + mR*R_t - aND*D_t*V_t - aVD*V_t*D_t + f*M,
            V_t = - mV*V_t -((aVR*D_t*R_t)/(1 + bR*R_t)) + aVD*V_t*D_t,
            P_t = - mP*P_t + ((aNP*N_t*P_t)/(1 + bP*P_t)) - aPH*P_t*H_t,
            H_t = - mH*H_t - ((aHR*H_t*R_t)/(1 + bR*R_t)) + aPH*P_t*H_t,
            R_t = - mR*R_t + ((aHR*H_t*R_t)/(1 + bR*R_t)) + ((aVR*D_t*R_t)/(1 + bR*R_t)) +
                (1 - f)*M)
        }
        )
    )

    return(list(output))
}


parms <- list(aND = 2.27e-7, aVD = 4.33e-7, aNP = 1e-3 * 6.12e-4,
        aPH = 1e-2 * 1.94e-5, aVR = 1e-4 * 2.15e-3, aHR = 1.13e-7,
        mN = 100/N0, mV = 2/V0, mP = 0.1/P0, mH = 0.001 * 2/H0, mR = 0.1 * 2/R0,
        G = 10, bP = 100, bR = 0,
        M_func = function(t) {ifelse(between(t,100,200),5,0)}, f = 2/3)

y <- c(N_t = 343000, D_t = 114000, V_t = 81, P_t = 4300, H_t = 24, R_t = 13)


output <- ode(y, 1:5000, diff_eq, parms)


output %>%
    as.data.frame %>%
    as_tibble %>%
    gather('pool', 'biomass', -time) %>%
    ggplot(aes(time, log(biomass), color = pool)) +
    geom_line(size = 1) +
    theme_classic() +
    scale_color_brewer(palette = 'Accent')
