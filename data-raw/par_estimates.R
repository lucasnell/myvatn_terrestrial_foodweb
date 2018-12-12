
library(readr)
library(dplyr)
library(purrr)

tibbles <- map(list.files("data-raw/parameters", "*.csv", full.names = TRUE),
               ~ read_csv(.x))
names(tibbles) <- gsub(".csv", "", list.files("data-raw/parameters", "*.csv"))

# Fix mistake in original file
tibbles$full_pars <- tibbles$full_pars %>% rename(mD = `\n            mD`)

par_estimates <- full_join(tibbles$full_pars,
                           tibbles$iN_eqs, by = "iN") %>%
    mutate(V = 1, H = 1, R = 1) %>%
    full_join(tibbles$structure_eqs %>% mutate(iN = 1000),
              by = c("V", "H", "R", "iN",
                     paste0(c("N", "D", "P", "V", "H", "R"), "eq"))) %>%
    select(V, H, R, iN, everything()) %>%
    arrange(desc(V), desc(H), desc(R), desc(iN)) %>%
    mutate_at(vars(matches("^l|^m|^a|^k")), mean, na.rm = TRUE) %>%
    mutate_at(vars(V, H, R), as.integer)

# usethis::use_data(par_estimates, overwrite = TRUE)
