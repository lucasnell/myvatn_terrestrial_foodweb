# =================================================================
# =================================================================

#       Create a pulse in midges through time.

# =================================================================
# =================================================================

#' Compute a midge pulse value for a given time point.
#'
#'
#' @noRd
#'
midge_pulse <- function(t, a, b, r, w, d) {
    max_adj <- sin(2 * pi * (0.75 - 1 / (3 * r)))
    max <- b * (1 + exp(-a * (1 + max_adj)))
    sines <- sin(2 * pi * (8 * t - 8 * d * w + 9 * r * w)/(12 * r * w))
    f <- max / (1+exp(a * (sines - max_adj)))
    return(f)
}


#' @describeIn food_web Test a midge pulse time series.
#'
#' @inheritParams food_web
#'
#' @return A numeric vector
#'
#' @export
#'
test_midges <- function(tmax, a, b, r, w, d, tstep = 1) {
    time <- seq(0, tmax, tstep)
    if (time[length(time)] < tmax) time <- c(time, tmax)
    return(midge_pulse(time, a, b, r, w, d))
}




# =================================================================
# =================================================================

#       Summarize a pulse in nitrogen content through time.

# =================================================================
# =================================================================



#' Summarize aspects of the pulse of nitrogen content in a time series.
#'
#' @param x A vector of nitrogen content through time.
#' @param measure Character specifying which summary measure to return.
#'     Options are as follows:
#'     \describe{
#'         \item{aamp}{absolute amplitude of the pulse.}
#'         \item{ramp}{relative amplitude of the pulse.}
#'         \item{peakt}{when the pulse peak occurs.}
#'         \item{len}{how long the pulse lasts.}
#'         \item{skew}{ratio of time before to after the pulse peak.}
#'     }
#' @param t0 Number of time points to subtract from the time to pulse peak.
#'     For example, if you want to ignore the first 10 time points, this would be
#'     set to `10`. This only affects the output if `aamp == "peakt"`.
#'     Defaults to `0`.
#' @param p Proportion of the total amplitude that is required for a time point to be
#'     considered inside a pulse.
#'     Only affects output (and is required) if `aamp %in% c("len", "skew")`.
#'     Defaults to `NULL`, which will return an error if the above condition is
#'     `TRUE`.
#'
#' @export
#'
pulse_summary <- function(x, measure, t0 = 0, p = NULL) {

    measure <- match.arg(measure, c("aamp", "ramp", "peakt", "len", "skew"))

    if (aamp == "aamp") {
        out <- pulse_amplitude(x, FALSE)
    } else if (aamp == "ramp") {
        out <- pulse_amplitude(x, TRUE)
    } else if (aamp == "peakt") {
        out <- pulse_peak_time(x, t0)
    } else if (aamp == "len") {
        if (is.null(p)) stop("\nspecify a p value for measure = \"len\"")
        out <- pulse_length(x, p)
    } else {
        if (is.null(p)) stop("\nspecify a p value for measure = \"skew\"")
        out <- pulse_skew(x, p)
    }
}


#' Compute the amplitude of the pulse.
#'
#' @inheritParams pulse_summary x
#' @param relative Logical for whether to return the relative amplitude (rather than
#'     absolute).
#'
#' @noRd
#'
pulse_amplitude <- function(x, relative) {
    z <- (max(x) - x[1])
    if (relative) z <- z / x[1]
    return(z)
}


#' Compute when the pulse peak occurs.
#'
#' @inheritParams pulse_summary x
#' @inheritParams pulse_summary t0
#'
#' @noRd
#'
pulse_peak_time <- function(x, t0) which(x == max(x)) - t0


#' Compute how long the pulse lasts.
#'
#' @inheritParams pulse_summary x
#' @inheritParams pulse_summary p
#'
#' @noRd
#'
pulse_length <- function(x, p) {
    amp <- max(x) - x[1]
    over_thresh <- which(x >= (p * amp + x[1]))
    if (length(over_thresh) == 0) {
        return(over_thresh)
    } else {
        return(diff(range(over_thresh)))
    }
}


#' Compute ratio of time before to after the pulse peak.
#'
#' @inheritParams pulse_summary x
#' @inheritParams pulse_summary p
#'
#' @noRd
#'
pulse_skew <- function(x, p) {
    max_t <- which(x == max(x))
    amp <- max(x) - x[1]
    over_thresh <- which(x >= (p * amp + x[1]))
    if (length(over_thresh) <= 3) stop("\np is too high, silly")
    before_peak <- max_t - over_thresh[1]
    after_peak <- tail(over_thresh, 1) - max_t
    return(before_peak / after_peak)
}
