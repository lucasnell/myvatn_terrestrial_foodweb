
#' Parameter estimates for the food web model.
#'
#' These estimates are computed separately when one or more pool is removed.
#'
#'
#' @format A data frame with 7 rows and 34 variables:
#' \describe{
#'   \item{Neq}{  Desired equilibrium value for N. }
#'   \item{Deq}{  Desired equilibrium value for D. }
#'   \item{Peq}{  Desired equilibrium value for P. }
#'   \item{Veq}{  Desired equilibrium value for V. }
#'   \item{Heq}{  Desired equilibrium value for H. }
#'   \item{Req}{  Desired equilibrium value for R. }
#'   \item{iN}{ Input to N. }
#'   \item{lD}{ Loss rates systems for D. }
#'   \item{lP}{ Loss rates systems for P. }
#'   \item{lV}{ Loss rates systems for V. }
#'   \item{lH}{ Loss rates systems for H. }
#'   \item{lR}{ Loss rates systems for R. }
#'   \item{lM}{ Loss rates systems for M. }
#'   \item{mN}{ Loss rates from pool N (returned to either N or D). }
#'   \item{mP}{ Loss rates from pool P (returned to either N or D). }
#'   \item{mD}{ Loss rates from pool D (returned to either N or D). }
#'   \item{mV}{ Loss rates from pool V (returned to either N or D). }
#'   \item{mH}{ Loss rates from pool H (returned to either N or D). }
#'   \item{mR}{ Loss rates from pool R (returned to either N or D). }
#'   \item{mM}{ Loss rates from pool M (returned to either N or D). }
#'   \item{kP}{ Carrying capacity for P. }
#'   \item{kV}{ Carrying capacity for V. }
#'   \item{kH}{ Carrying capacity for H. }
#'   \item{kR}{ Carrying capacity for R. }
#'   \item{aNP}{ Uptake rate for NP. }
#'   \item{aDV}{ Uptake rate for DV. }
#'   \item{aPH}{ Uptake rate for PH. }
#'   \item{aR}{ Uptake rate for R. }
#' }
#'
"par_estimates"
