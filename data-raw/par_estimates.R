
library(readr)
library(dplyr)

par_estimates <- read_csv("data-raw/par_estimates.csv",
                          col_types = do.call(cols, as.list(c(rep("i", 3),
                                                              rep("d", 31)))))

usethis::use_data(par_estimates)
