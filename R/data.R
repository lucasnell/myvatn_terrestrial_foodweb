
#' Parameter estimates for the food web model.
#'
#' These estimates are computed separately when one or more pool is removed.
#'
#'
#' @format A data frame with 8 rows and 37 variables:
#' \describe{
#'   \item{V}{ Integer for whether the `V` pool is included (0 = `FALSE`, 1 = `TRUE`). }
#'   \item{H}{ Integer for whether the `H` pool is included (0 = `FALSE`, 1 = `TRUE`). }
#'   \item{R}{ Integer for whether the `R` pool is included (0 = `FALSE`, 1 = `TRUE`). }
#'   \item{iN}{ Input to N. }
#'   \item{lD}{ Loss rates systems for D. }
#'   \item{lP}{ Loss rates systems for P. }
#'   \item{lV}{ Loss rates systems for V. }
#'   \item{lH}{ Loss rates systems for H. }
#'   \item{lR}{ Loss rates systems for R. }
#'   \item{lM}{ Loss rates systems for M. }
#'   \item{mN}{ Loss rates from pool N (returned to either N or D). }
#'   \item{mD}{ Loss rates from pool D (returned to either N or D). }
#'   \item{mP0}{ Baseline Loss rates from pool P (returned to either N or D). }
#'   \item{mD0}{ Baseline Loss rates from pool D (returned to either N or D). }
#'   \item{mV0}{ Baseline Loss rates from pool V (returned to either N or D). }
#'   \item{mH0}{ Baseline Loss rates from pool H (returned to either N or D). }
#'   \item{mR0}{ Baseline Loss rates from pool R (returned to either N or D). }
#'   \item{mP}{ DD Loss rates from pool P (returned to either N or D). }
#'   \item{mV}{ DD Loss rates from pool V (returned to either N or D). }
#'   \item{mH}{ DD Loss rates from pool H (returned to either N or D). }
#'   \item{mR}{ DD Loss rates from pool R (returned to either N or D). }
#'   \item{mM}{ Loss rates from pool M (returned to either N or D). }
#'   \item{aNP}{ Uptake rate for NP. }
#'   \item{aDV}{ Uptake rate for DV. }
#'   \item{aPH}{ Uptake rate for PH. }
#'   \item{aR}{ Uptake rate for R. }
#'   \item{aM}{ Uptake rate of M for R. }
#'   \item{hN}{ Handling time for uptake of N. }
#'   \item{hD}{ Handling time for uptake of D. }
#'   \item{hP}{ Handling time for uptake of P. }
#'   \item{hVHM}{ Handling time for uptake of V, H, and M. }
#'   \item{Neq}{  Equilibrium value for N. }
#'   \item{Deq}{  Equilibrium value for D. }
#'   \item{Peq}{  Equilibrium value for P. }
#'   \item{Veq}{  Equilibrium value for V. }
#'   \item{Heq}{  Equilibrium value for H. }
#'   \item{Req}{  Equilibrium value for R. }
#' }
#'
"par_estimates"
