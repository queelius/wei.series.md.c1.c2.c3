#' Likelihood model for Weibull series systems from masked data.
#'
#' Functions include the log-likelihood, score, and hessian of the log-likelihood
#' functions. We provide analytical solutions to the log-likelihood and score
#' functions. The hessian is computed numerically by taking the Jacobian of the
#' score function using the `jacobian` function in the R package `numDeriv`.
#' 
#' The likelihood model is given by two types of data:
#' 1. Masked component cause of failure data with exact failure time
#' 2. Right-censored system lifetime data
#' 
#' Masked component data approximately satisfies the following conditions:
#' C1: Pr(K in C) = 1
#' C2: Pr(C=c | K=j, T=t) = Pr(C=c | K=j', T=t)
#'     for any j, j' in c.
#' C3: masking probabilities are independent of theta
#' 
#' @author Alex Towell
#' @name Weibull series MLE
#' @keywords weibull, distribution, series, statistics, masked data
NULL

#' Generates a log-likelihood function for a Weibull series system with respect
#' to parameter `theta` (shape, scale) for masked data with candidate sets
#' that satisfy conditions C1, C2, and C3 and right-censored data.
#'
#' @param df (masked) data frame
#' @param theta parameter vector (shape1, scale1, ..., shapem, scalem)
#' @param candset prefix of Boolean matrix encoding of candidate sets,
#'       defaults to `x`, e.g., `x1,...,xm`.
#' @param lifetime system lifetime (optionally right-censored) column name
#' @param right_censoring_indicator right-censoring indicator column name, if
#'       TRUE, then the system lifetime is observed, otherwise it is right-censored.
#'      If there is no right-censoring indicator column by the given name, then the
#'     system lifetimes are assumed to be observed.
#' @returns Log-likelihood with respect to `theta` given `df`
#' @importFrom md.tools md_decode_matrix
#' @export
loglik_wei_series_md_c1_c2_c3 <- function(
    df,
    theta,
    candset = "x",
    lifetime = "t",
    right_censoring_indicator = "delta") {
     
    if (!lifetime %in% colnames(df)) {
        stop("lifetime variable not in colnames(df)")
    }

    n <- nrow(df)
    if (n < 1) {
        stop("sample size must be greater than 0")
    }

    C <- md_decode_matrix(df, candset)
    if (is.null(C)) {
        stop("no candidate set found for candset")
    }
    m <- ncol(C)

    if (right_censoring_indicator %in% colnames(df)) {
        delta <- df[[right_censoring_indicator]]
    } else {
        delta <- rep(TRUE, n)
    }

    t <- df[[lifetime]]

    k <- length(theta)
    stopifnot(k == 2 * m)
    shapes <- theta[seq(1, k, 2)]
    scales <- theta[seq(2, k, 2)]

    s <- 0
    for (i in 1:n) {
        s <- s - sum((t[i] / scales)^shapes)
        if (delta[i]) {
            s <- s + log(sum(shapes[C[i, ]] / scales[C[i, ]] *
                (t[i] / scales[C[i, ]])^(shapes[C[i, ]] - 1)))
        }
    }
    s
}

#' Slowly computes the score function (gradient of the log-likelihood function) for a
#' Weibull series system with respect to parameter `theta` (shape, scale) for masked
#' data with candidate sets that satisfy conditions C1, C2, and C3 and right-censored
#' data. This is simple loop-based implementation of the score function, which is
#' easier to verify the correctness of.
#'
#' @param df (masked) data frame
#' @param theta parameter vector (shape1, scale1, ..., shapem, scalem)
#' @param candset prefix of Boolean matrix encoding of candidate sets,
#'        defaults to `x`, e.g., `x1,...,xm`.
#' @param lifetime system lifetime (optionally right-censored) column name
#' @param right_censoring_indicator right-censoring indicator column name, if
#'        TRUE, then the system lifetime is observed, otherwise it is right-censored.
#'        If there is no right-censoring indicator column by the given name, then the
#'        system lifetimes are assumed to be observed.
#' @returns Score with respect to `theta` given `df`
#' @export
score_wei_series_md_c1_c2_c3_slow <- function(
    df,
    theta,
    candset = "x",
    lifetime = "t",
    right_censoring_indicator = "delta") {

    if (!lifetime %in% colnames(df)) {
        stop("lifetime variable not in colnames(df)")
    }

    n <- nrow(df)
    if (n < 1) {
        stop("sample size must be greater than 0")
    }

    C <- md_decode_matrix(df, candset)
    if (is.null(C)) {
        stop("no candidate set found for candset")
    }
    m <- ncol(C)

    if (right_censoring_indicator %in% colnames(df)) {
        delta <- df[[right_censoring_indicator]]
    } else {
        delta <- rep(TRUE, n)
    }

    t <- df[[lifetime]]
    shapes <- theta[seq(1, length(theta), 2)]
    scales <- theta[seq(2, length(theta), 2)]
    shape_scores <- rep(0, m)
    scale_scores <- rep(0, m)

    for (j in 1:m) {
        # let's do shape scores
        for (i in 1:n) {
            rt.term <- -(t[i] / scales[j])^shapes[j] * log(t[i] / scales[j])
            mask.term <- 0
            if (delta[i] && C[i,j]) {
                denom <- 0
                for (k in 1:m) {
                    if (C[i,k]) {
                        denom <- denom + shapes[k] / scales[k] *
                            (t[i] / scales[k])^(shapes[k] - 1)
                    }
                }
                numer <- 1/t[i] * (t[i] / scales[j])^shapes[j] *
                    (1 + shapes[j] * log(t[i] / scales[j]))
                mask.term <- numer / denom
            }
            shape_scores[j] <- shape_scores[j] + rt.term + mask.term
        }

        # let's do scale scores
        for (i in 1:n) {
            rt.term <- (shapes[j] / scales[j]) * (t[i] / scales[j])^shapes[j]
            mask.term <- 0
            if (delta[i] && C[i,j]) {
                denom <- 0
                for (k in 1:m) {
                    if (C[i, k]) {
                        denom <- denom + shapes[k] / scales[k] *
                            (t[i] / scales[k])^(shapes[k] - 1)
                    }
                }

                numer <- (shapes[j] / scales[j])^2 * (t[i] / scales[j])^(shapes[j] - 1)
                mask.term <- numer / denom
            }
            scale_scores[j] <- scale_scores[j] + rt.term - mask.term
        }
    }
    
    scr <- rep(0, length(theta))
    scr[seq(1, length(theta), 2)] <- shape_scores
    scr[seq(2, length(theta), 2)] <- scale_scores
    scr
}

#' Generates a hessian of the log-likelihood function (negative of the observed
#' FIM) for a Weibull series system with respect to parameter `theta` (shape, scale)
#' for masked data with candidate sets that satisfy conditions C1, C2, and C3 and
#' right-censored data.
#'
#' @param df (masked) data frame
#' @param theta parameter vector (shape1, scale1, ..., shapem, scalem)
#' @param candset prefix of Boolean matrix encoding of candidate sets,
#'        defaults to `x`, e.g., `x1,...,xm`.
#' @param lifetime system lifetime (optionally right-censored) column name
#' @param right_censoring_indicator right-censoring indicator column name, if
#'        TRUE, then the system lifetime is observed, otherwise it is right-censored.
#'        If there is no right-censoring indicator column by the given name, then the
#'        system lifetimes are assumed to be observed.
#' @param ... additional arguments passed to `method.args` in `hessian`
#' @returns Score with respect to `theta` given `df`
#' @importFrom numDeriv jacobian
#' @export
hessian_wei_series_md_c1_c2_c3 <- function(
    df,
    theta,
    candset = "x",
    lifetime = "t",
    right_censoring_indicator = "delta",
    ...) {

    jacobian(
        func = score_wei_series_md_c1_c2_c3,
        x = theta,
        df = df,
        candset = candset,
        lifetime = lifetime,
        right_censoring_indicator = right_censoring_indicator,
        method.args = list(r = 6, ...))
}

#' Cmputes the score function (gradient of the log-likelihood function) for a
#' Weibull series system with respect to parameter `theta` (shape, scale) for masked
#' data with candidate sets that satisfy conditions C1, C2, and C3 and right-censored
#' data.
#'
#' @param df (masked) data frame
#' @param theta parameter vector (shape1, scale1, ..., shapem, scalem)
#' @param candset prefix of Boolean matrix encoding of candidate sets,
#'        defaults to `x`, e.g., `x1,...,xm`.
#' @param lifetime system lifetime (optionally right-censored) column name
#' @param right_censoring_indicator right-censoring indicator column name, if
#'        TRUE, then the system lifetime is observed, otherwise it is right-censored.
#'        If there is no right-censoring indicator column by the given name, then the
#'        system lifetimes are assumed to be observed.
#' @param ... additional arguments passed to `method.args` in `grad`
#' @returns Score with respect to `theta` given `df`
#' @export
score_wei_series_md_c1_c2_c3 <- function(
    df,
    theta,
    candset = "x",
    lifetime = "t",
    right_censoring_indicator = "delta",
    ...) {

    if (!lifetime %in% colnames(df)) {
        stop("lifetime variable not in colnames(df)")
    }

    n <- nrow(df)
    if (n < 1) {
        stop("sample size must be greater than 0")
    }

    C <- md_decode_matrix(df, candset)
    if (is.null(C)) {
        stop("no candidate set found for candset")
    }
    m <- ncol(C)

    if (right_censoring_indicator %in% colnames(df)) {
        delta <- df[[right_censoring_indicator]]
    } else {
        delta <- rep(TRUE, n)
    }

    t <- df[[lifetime]]
    shapes <- theta[seq(1, length(theta), 2)]
    scales <- theta[seq(2, length(theta), 2)]
    shape_scores <- rep(0, m)
    scale_scores <- rep(0, m)

    for (i in 1:n) {
        rt.term.shapes <- -(t[i] / scales)^shapes * log(t[i] / scales)
        rt.term.scales <- (shapes / scales) * (t[i] / scales)^shapes

        # Initialize mask terms to 0
        mask.term.shapes <- rep(0, m)
        mask.term.scales <- rep(0, m)
        
        # Perform the calculations only where C[i,] is TRUE
        if (delta[i]) {
            cindex <- C[i,]
            denom <- sum(shapes[cindex] / scales[cindex] * (t[i] / scales[cindex])^(shapes[cindex] - 1))
            
            numer.shapes <- 1/t[i] * (t[i] / scales[cindex])^shapes[cindex] *
                (1 + shapes[cindex] * log(t[i] / scales[cindex]))
            mask.term.shapes[cindex] <- numer.shapes / denom

            numer.scales <- (shapes[cindex] / scales[cindex])^2 * (t[i] / scales[cindex])^(shapes[cindex] - 1)
            mask.term.scales[cindex] <- numer.scales / denom
        }

        shape_scores <- shape_scores + rt.term.shapes + mask.term.shapes
        scale_scores <- scale_scores + rt.term.scales - mask.term.scales
    }

    scr <- rep(0, length(theta))
    scr[seq(1, length(theta), 2)] <- shape_scores
    scr[seq(2, length(theta), 2)] <- scale_scores
    scr
}

