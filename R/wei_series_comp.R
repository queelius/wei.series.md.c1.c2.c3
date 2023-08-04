#' Weibull Series System Hazard Function
#' @param t time
#' @param shape shape parameter
#' @param scale scale parameter
#' @export
hazard_wei <- Vectorize(function(t, shape, scale) {
    shape / scale * (t / scale)^(shape - 1)
}, vectorize.args = "t")

#' Mean-Time-To-Failure for Weibull
#' @param shape shape parameter
#' @param scale scale parameter
#' @export
wei_mttf <- function(shape, scale) {
    scale * gamma(1 + 1 / shape)
}
