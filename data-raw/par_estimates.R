
library(readr)
library(dplyr)
library(purrr)

tibbles <- map(list.files("data-raw/parameters", "*.csv", full.names = TRUE),
               ~ read_csv(.x))
names(tibbles) <- gsub(".csv", "", list.files("data-raw/parameters", "*.csv"))

# Function to fill NAs in columns with only 1 non-NA value.
# Other columns are unchanged
fill_nas <- function(x) {
    if (sum(is.na(x)) == length(x) - 1) {
        return(rep(x[!is.na(x)], length(x)))
    }
    return(x)
}

par_estimates <- full_join(tibbles$full_pars,
                           tibbles$iN_eqs, by = "iN") %>%
    mutate(V = 1, H = 1, R = 1) %>%
    full_join(tibbles$structure_eqs %>% mutate(iN = 10),
              by = c("V", "H", "R", "iN",
                     paste0(c("N", "D", "P", "V", "H", "R"), "eq"))) %>%
    select(V, H, R, iN, everything()) %>%
    arrange(desc(V), desc(H), desc(R), desc(iN)) %>%
    mutate_all(fill_nas) %>%
    mutate_at(vars(V, H, R), as.integer)


#'
#' Below changes attack rate on midges to the relative attack rate in relation to
#' the rate for herbivores and detritivores.
#'
#' The overall attack rate was used to originally calculate the equilibrium values,
#' but we've decided to use the relative rate in the paper because it improves
#' the clarity of discussion.
#'
par_estimates <- par_estimates %>%
    mutate(aM = aM / aR) %>%
    rename(q = aM)


sw <- function(x) starts_with(x, ignore.case = FALSE)
ew <- function(x) ends_with(x, ignore.case = FALSE)

#'
#' Below allows midges to have their own handling time.
#'
par_estimates <- par_estimates %>%
    mutate(hM = hVHM) %>%
    rename(hVH = hVHM) %>%
    # To adjust for new parameter names in manuscript:
    rename(iI = iN,
           aIP = aNP,
           hI = hN,
           X = R,
           lX = lR,
           mX = mR,
           aX = aR,
           hX = hVH,
           muI = mN,
           muD = mD,
           muP = mP0,
           muV = mV0,
           muH = mH0,
           muX = mR0,
           muM = mM,
           Ieq = Neq,
           Xeq = Req) %>%
    select(V, H, X,
           iI,
           lD, lP, lV, lH, lX, lM,
           muI, muD, muP, muV, muH, muX, muM,
           mP, mV, mH, mX,
           q,
           aIP, aDV, aPH, aX,
           hI, hD, hP, hX, hM,
           Ieq, Deq, Peq, Veq, Heq, Xeq)





# usethis::use_data(par_estimates, overwrite = TRUE)
