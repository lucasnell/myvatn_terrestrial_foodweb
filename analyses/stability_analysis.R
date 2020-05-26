##### Load packages and data

library(mtf)
library(tidyverse)

pars <- pars <- par_estimates %>%
    filter(V == 1, X == 1, H == 1, iI == 10) %>%
    dplyr::select(-V, -X, -H)
pars <- as.list(pars)





##### Define functions for derivatives

dI_fn <- deriv(~iI - aIP*I*P/(1 + aIP*hI*I) + (1 - lD)*muD*D - muI*I,
            c("I", "D", "P", "V", "H", "X"),
            function(I, D, P, V, H, X, M,
                     iI, lD, lP, lV, lH, lX, lM, muI, muD, muP, muV, muH, muX, muM, mP, mV,
                     mH, mX, q,  aIP, aDV, aPH, aX, hI, hD, hP, hX, hM){})

dD_fn <- deriv(~(1 - lP)*(muP + mP*P)*P + (1 - lV)*(muV + mV*V)*V + (1 - lH)*(muH + mH*H)*H +
                (1 - lX)*(muX + mX*X)*X + (1 - lM)*muM*M - aDV*D*V/(1 + aDV*hD*D) - muD*D,
            c("I", "D", "P", "V", "H", "X"),
            function(I, D, P, V, H, X, M,
                     iI, lD, lP, lV, lH, lX, lM, muI, muD, muP, muV, muH, muX, muM, mP, mV,
                     mH, mX, q,  aIP, aDV, aPH, aX, hI, hD, hP, hX, hM){})

dP_fn <- deriv(~aIP*I*P/(1 + aIP*hI*I) - aPH*P*H/(1 + aPH*hP*P) - (muP + mP*P)*P,
            c("I", "D", "P", "V", "H", "X"),
            function(I, D, P, V, H, X, M,
                     iI, lD, lP, lV, lH, lX, lM, muI, muD, muP, muV, muH, muX, muM, mP, mV,
                     mH, mX, q,  aIP, aDV, aPH, aX, hI, hD, hP, hX, hM){})

dV_fn <- deriv(~aDV*D*V/(1 + aDV*hD*D) - (aX*V*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M) -
                (muV + mV*V)*V,
            c("I", "D", "P", "V", "H", "X"),
            function(I, D, P, V, H, X, M,
                     iI, lD, lP, lV, lH, lX, lM, muI, muD, muP, muV, muH, muX, muM, mP, mV,
                     mH, mX, q,  aIP, aDV, aPH, aX, hI, hD, hP, hX, hM){})

dH_fn <- deriv(~aPH*P*H/(1 + aPH*hP*P) - (aX*H*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M) -
                (muH + mH*H)*H,
            c("I", "D", "P", "V", "H", "X"),
            function(I, D, P, V, H, X, M,
                     iI, lD, lP, lV, lH, lX, lM, muI, muD, muP, muV, muH, muX, muM, mP, mV,
                     mH, mX, q,  aIP, aDV, aPH, aX, hI, hD, hP, hX, hM){})

dX_fn <- deriv(~(aX*V*X + aX*H*X + (aX * q)*M*X)/(1 + aX*hX*(V + H) + (aX * q)*hM*M) -
                (muX + mX*X)*X,
            c("I", "D", "P", "V", "H", "X"),
            function(I, D, P, V, H, X, M,
                     iI, lD, lP, lV, lH, lX, lM, muI, muD, muP, muV, muH, muX, muM, mP, mV,
                     mH, mX, q,  aIP, aDV, aPH, aX, hI, hD, hP, hX, hM){})




##### Evaluate derivatives at equilibrium

dI <- with(pars,
        dI_fn(Ieq, Deq, Peq, Veq, Heq, Xeq, 0,
            iI, lD, lP, lV, lH, lX, lM, muI, muD, muP, muV, muH, muX, muM, mP, mV,
            mH, mX, q,  aIP, aDV, aPH, aX, hI, hD, hP, hX, hM))

dD <- with(pars,
        dD_fn(Ieq, Deq, Peq, Veq, Heq, Xeq, 0,
           iI, lD, lP, lV, lH, lX, lM, muI, muD, muP, muV, muH, muX, muM, mP, mV,
           mH, mX, q,  aIP, aDV, aPH, aX, hI, hD, hP, hX, hM))

dP <- with(pars,
        dP_fn(Ieq, Deq, Peq, Veq, Heq, Xeq, 0,
           iI, lD, lP, lV, lH, lX, lM, muI, muD, muP, muV, muH, muX, muM, mP, mV,
           mH, mX, q,  aIP, aDV, aPH, aX, hI, hD, hP, hX, hM))

dV <- with(pars,
        dV_fn(Ieq, Deq, Peq, Veq, Heq, Xeq, 0,
              iI, lD, lP, lV, lH, lX, lM, muI, muD, muP, muV, muH, muX, muM, mP, mV,
              mH, mX, q,  aIP, aDV, aPH, aX, hI, hD, hP, hX, hM))

dH <- with(pars,
        dH_fn(Ieq, Deq, Peq, Veq, Heq, Xeq, 0,
              iI, lD, lP, lV, lH, lX, lM, muI, muD, muP, muV, muH, muX, muM, mP, mV,
              mH, mX, q,  aIP, aDV, aPH, aX, hI, hD, hP, hX, hM))

dX <- with(pars,
        dX_fn(Ieq, Deq, Peq, Veq, Heq, Xeq, 0,
              iI, lD, lP, lV, lH, lX, lM, muI, muD, muP, muV, muH, muX, muM, mP, mV,
              mH, mX, q,  aIP, aDV, aPH, aX, hI, hD, hP, hX, hM))




##### Define jacobian

jac <- t(sapply(list(dI, dD, dP, dV, dH, dX),
              function(x){attributes(x)$gradient}))
colnames(jac) <- c("I", "D", "P", "V", "H", "X")
rownames(jac) <- c("I", "D", "P", "V", "H", "X")





##### Eigenvalues

eigen_analysis <- eigen(jac)
eigen_values <- eigen_analysis$values
leading_eign <- max(Re(eigen_values))

leading_eign # leading eigenvalue is negative, so equilibrium is stable
eigen_values # imaginary eigenvalues implies some kind of cyclic approach to equilibrium
# (but this gets complicated)
