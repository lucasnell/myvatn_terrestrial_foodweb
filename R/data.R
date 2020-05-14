
#' Parameter estimates for the food web model.
#'
#' These estimates are computed separately when one or more pool is removed.
#'
#'
#' @format A data frame with 8 rows and 37 variables:
#' \describe{
#'   \item{V}{ Integer for whether the `V` pool is included (0 = `FALSE`, 1 = `TRUE`). }
#'   \item{H}{ Integer for whether the `H` pool is included (0 = `FALSE`, 1 = `TRUE`). }
#'   \item{X}{ Integer for whether the `X` pool is included (0 = `FALSE`, 1 = `TRUE`). }
#'   \item{iI}{ Input to I. }
#'   \item{lD}{ Loss rates systems for D. }
#'   \item{lP}{ Loss rates systems for P. }
#'   \item{lV}{ Loss rates systems for V. }
#'   \item{lH}{ Loss rates systems for H. }
#'   \item{lX}{ Loss rates systems for X. }
#'   \item{lM}{ Loss rates systems for M. }
#'   \item{muI}{ Density-independent loss rate from pool I (returned to either I or D). }
#'   \item{muD}{ Density-independent loss rate from pool D (returned to either I or D). }
#'   \item{muP}{ Density-independent loss rate from pool P (returned to either I or D). }
#'   \item{muV}{ Density-independent loss rate from pool V (returned to either I or D). }
#'   \item{muH}{ Density-independent loss rate from pool H (returned to either I or D). }
#'   \item{muX}{ Density-independent loss rate from pool X (returned to either I or D). }
#'   \item{muM}{ Density-independent loss rate from pool M (returned to either I or D). }
#'   \item{mP}{ Density-dependent loss rate from pool P (returned to either I or D). }
#'   \item{mV}{ Density-dependent loss rate from pool V (returned to either I or D). }
#'   \item{mH}{ Density-dependent loss rate from pool H (returned to either I or D). }
#'   \item{mX}{ Density-dependent loss rate from pool X (returned to either I or D). }
#'   \item{q}{ Relative rate of predator attack on midges in relation to that on herbivores and detritivores. }
#'   \item{aIP}{ Uptake rate for IP. }
#'   \item{aDV}{ Uptake rate for DV. }
#'   \item{aPH}{ Uptake rate for PH. }
#'   \item{aX}{ Uptake rate for X. }
#'   \item{hI}{ Handling time for uptake of I. }
#'   \item{hD}{ Handling time for uptake of D. }
#'   \item{hP}{ Handling time for uptake of P. }
#'   \item{hX}{ Handling time for uptake of V and H by X. }
#'   \item{hM}{ Handling time for uptake of M. }
#'   \item{Ieq}{  Equilibrium value for I. }
#'   \item{Deq}{  Equilibrium value for D. }
#'   \item{Peq}{  Equilibrium value for P. }
#'   \item{Veq}{  Equilibrium value for V. }
#'   \item{Heq}{  Equilibrium value for H. }
#'   \item{Xeq}{  Equilibrium value for X. }
#' }
#'
"par_estimates"



