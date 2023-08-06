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
#' @param shapes vector of weibull shape parameters for weibull lifetime
#'               components
#' @param scales vector of weibull scale parameters for weibull lifetime
#'               components
#' @param tol stopping condition, default is 1e-5
#' @param t0 initial guess, default is 1
#' @export
qwei_series <- Vectorize(function(p, shapes, scales, tol=1e-5, t0 = 1) {
    stopifnot(length(scales) == length(shapes),
        all(shapes > 0), all(scales > 0))

    if (p <= 0) {
        return(0)
    }
    if (p >= 1) {
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
        if (abs(t1 - t0) < tol) {
            break
        }
        t0 <- t1
    }
    return(t1)
}, vectorize.args = "p")

#' Sampler for weibull series.
#'
#' @param n sample size
#' @param shapes shape parameters for weibull component lifetimes
#' @param scales scale parameters for weibull component lifetimes
#' @importFrom stats rweibull
#' @export
rwei_series <- function(n, shapes, scales) {
    stopifnot(n > 0)
    m <- length(scales)
    stopifnot(m == length(shapes), all(shapes > 0), all(scales > 0))

    t <- matrix(nrow = n, ncol = m)
    for (j in 1:m)
        t[, j] <- rweibull(n, shape = shapes[j], scale = scales[j])
    apply(t, 1, min)
}

#' pdf for weibull series
#'
#' @param t series system lifetime
#' @param shapes shape parameters for weibull component lifetimes
#' @param scales scale parameters for weibull component lifetimes
#' @export
dwei_series <- Vectorize(function(t, shapes, scales) {
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
#' @param shapes shape parameters for weibull component lifetimes
#' @param scales scale parameters for weibull component lifetimes
#' @export
surv_wei_series <- Vectorize(function(t, shapes, scales) {
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
#' @param shapes shape parameters for weibull component lifetimes
#' @param scales scale parameters for weibull component lifetimes
#' @export
hazard_wei_series <- Vectorize(function(t, shapes, scales) {
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))
    
    ifelse(t <= 0, NA, sum(shapes/scales*(t/scales)^(shapes-1)))
}, vectorize.args="t")

#' The cumulative distribution function for Weibull series
#'
#' @param t series system lifetime
#' @param scales scale parameters for Weibull component lifetimes
#' @param shapes shape parameters for Weibull component lifetimes
#' @export
pwei_series <- Vectorize(function(t, shapes, scales) {
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))
    ifelse(t < 0, 0, 1-exp(-sum((t/scales)^shapes)))
}, vectorize.args = "t")

#' The mean time to failure for the Weibull series system.
#' 
#' @param shapes shape parameters for Weibull component lifetimes
#' @param scales scale parameters for Weibull component lifetimes
#' @importFrom stats integrate
#' @export
wei_series_mttf <- function(shapes, scales) {
    stopifnot(length(shapes) == length(scales),
        all(shapes > 0), all(scales > 0))
    integrate(surv_wei_series, shapes = shapes, scales = scales, 0, Inf)$value
}

#' The component cause of failure for the Weibull series distribution.
#' 
#' @param k component index
#' @param t lifetime of system, NA if unconditional distribution of K[i]
#' @param shapes shape parameters for Weibull component lifetimes
#' @param scales scale parameters for Weibull component lifetimes
#' @param mc logical, if TRUE, use monte carlo integration instead, which
#'           may be more robust for some parameters, particularly small
#'           shape parameters
#' @return the probability that component k is the cause of failure,
#'         for the unconditional random variable K[i]
#' @export
wei_series_cause <- Vectorize(function(k, shapes, scales, t = NA, mc = FALSE, n = 10000L) {
    m <- length(shapes)
    k <- as.integer(k)
    stopifnot(m == length(scales))
    if (k < 1 || k > m) {
        return(0)
    }

    if (!is.na(t)) {
        if (t < 0) {
            return(NA)
        }
        return(hazard_wei(t = t, shape = shapes[k], scale = scales[k]) /
            hazard_wei_series(t = t, shapes = shapes, scales = scales))
    } else {
        # unconditional, E(h_j(T) / h(T))
        if (mc) {
            ts <- rwei_series(n, shapes, scales)
            return(mean(x = 
                hazard_wei(t = ts, shape = shapes[k], scale = scales[k]) /
                hazard_wei_series(t = ts, shapes = shapes, scales = scales),
                na.rm = TRUE))
        } else {
            return(integrate(function(t) {
                dwei_series(t, shapes = shapes, scales = scales) *
                    hazard_wei(t = t, shape = shapes[k], scale = scales[k]) /
                    hazard_wei_series(t = t, shapes = shapes, scales = scales) },
                0, Inf)$value)
        }
    }
}, vectorize.args = "k")
