#' Weibull series system
#'
#' This file contains functions related to the Weibull series distribution.
#' Functions include simulation, pdf, cdf, quantile, and other related
#' functions for the Weibull series distribution.
#' 
#' @author Alex Towell
#' @name Weibull series
#' @keywords weibull, distribution, series, statistics
NULL

#' Quantile function (inverse of the cdf).
#' By definition, the quantile p * 100% for a strictly monotonically increasing
#' cdf F is the value t that satisfies \code{F(t) - p = 0}.
#' We solve for t using newton's method.
#'
#' @param p vector of probabilities.
#' @param scales vector of weibull scale parameters for weibull lifetime
#'               components
#' @param shapes vector of weibull shape parameters for weibull lifetime
#'               components
#' @param eps stopping condition, default is 1e-5
#' @param t0 initial guess, default is 1
#' @export
qwei_series <- Vectorize(function(p, scales, shapes, eps=1e-5, t0 = 1) {
    stopifnot(length(scales) == length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    if (p < 0 || p > 1) {
        stop("p must be in [0,1]")
    }
    if (p == 0) {
        return(0)
    }
    if (p == 1) {
        return(Inf)
    }

    t1 <- NULL
    repeat {
        alpha <- 1
        repeat {
            t1 <- t0 - alpha * (sum((t0/scales)^shapes) + log(1-p)) /
                sum(shapes*t0^(shapes-1)/scales^shapes)
            if (!is.nan(t1) && t1 > 0) {
                break
            }
            alpha <- alpha / 2
        }
        if (abs(t1 - t0) < eps) {
            break
        }
        t0 <- t1
    }
    t1
}, vectorize.args="p")

#' Sampler for weibull series.
#'
#' @param n sample size
#' @param scales scale parameters for weibull component lifetimes
#' @param shapes shape parameters for weibull component lifetimes
#' @importFrom stats rweibull
#' @export
rwei_series <- function(n,scales,shapes) {
    stopifnot(n > 0)
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    t <- matrix(nrow=n,ncol=m)
    for (j in 1:m)
        t[,j] <- rweibull(n,scale=scales[j],shape=shapes[j])
    apply(t,1,min)
}

#' pdf for weibull series
#'
#' @param t series system lifetime
#' @param scales scale parameters for weibull component lifetimes
#' @param shapes shape parameters for weibull component lifetimes
#' @export
dwei_series <- Vectorize(function(t,scales,shapes) {
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    d <- ifelse(t < 0,
                0,
                sum(shapes/scales*(t/scales)^(shapes-1))*exp(-sum((t/scales)^shapes)))
    ifelse(is.nan(d), 0, d)
}, vectorize.args="t")

#' Survival function for weibull series
#'
#' @param t series system lifetime
#' @param scales scale parameters for weibull component lifetimes
#' @param shapes shape parameters for weibull component lifetimes
#' @export
surv_wei_series <- Vectorize(function(t,scales,shapes) {
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    ifelse(t < 0,
           1,
           exp(-sum((t/scales)^shapes)))
}, vectorize.args="t")

#' Hazard function for weibull series.
#'
#' @param t series system lifetime
#' @param scales scale parameters for weibull component lifetimes
#' @param shapes shape parameters for weibull component lifetimes
#' @export
hazard_wei_series <- Vectorize(function(t,scales,shapes) {
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    ifelse(t < 0,
           0,
           sum(shapes/scales*(t/scales)^(shapes-1)))
}, vectorize.args="t")


#' The cumulative distribution function for Weibull series
#'
#' @param t series system lifetime
#' @param scales scale parameters for Weibull component lifetimes
#' @param shapes shape parameters for Weibull component lifetimes
#' @export
pwei_series <- Vectorize(function(t,scales,shapes) {
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))
    ifelse(t < 0, 0, 1-exp(-sum((t/scales)^shapes)))
}, vectorize.args="t")
