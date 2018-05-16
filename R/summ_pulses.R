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
#'     * `aamp`: absolute amplitude of the pulse.
#'     * `ramp`: relative amplitude of the pulse.
#'     * `peakt`: when the pulse peak occurs.
#'     * `len`: how long the pulse lasts.
#'     * `skew`: ratio of time before to after the pulse peak.
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
