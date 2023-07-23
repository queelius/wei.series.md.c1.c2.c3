#' Slow L-BFGS-B MLE solver for the Weibull series model with
#' masked component cause of failure with observation of exact system lifetime and
#' right censoring.
#' @param df data frame, right-censored lifetimes with masked component cause of failure
#' @param theta0 initial parameter vector
#' @param ... additional arguments passed to `loglik_wei_series_md_c1_c2_c3` and
#'            `score_wei_series_md_c1_c2_c3`.
#' @param control list of control parameters passed to `optim` for `method == "L-BFGS-B"`
#' @return list with components:
#' - `par` final parameter vector
#' - `value` final log-likelihood
#' - `counts` number of function evaluations, gradient evaluations, and Hessian evaluations
#' - `convergence` convergence code
#' - `message` convergence message
#' - `hessian` estimated Hessian
#' @importFrom stats optim
#' @export
mle_lbfgsb_wei_series_md_c1_c2_c3_slow <- function(
    df,
    theta0,
    hessian = TRUE,
    ...,
    control = list()) {

    defaults <- list(fnscale = -1)
    control <- modifyList(defaults, control)

    optim(
        par = theta0,
        fn = function(theta) {
            loglik_wei_series_md_c1_c2_c3(df = df, theta = theta, ...)
        },
        gr = function(theta) {
            score_wei_series_md_c1_c2_c3_slow(df = df, theta = theta, ...)
        },
        hessian = hessian,
        lower = rep(1e-3, length(theta0)),
        method = "L-BFGS-B",
        control = control)
}

#' L-BFGS-B MLE solver for the Weibull series model with
#' masked component cause of failure with observation of exact system lifetime and
#' right censoring.
#' @param df data frame, right-censored lifetimes with masked component cause of failure
#' @param theta0 initial parameter vector
#' @param ... additional arguments passed to `loglik_wei_series_md_c1_c2_c3` and
#'            `score_wei_series_md_c1_c2_c3`.
#' @param control list of control parameters passed to `optim` for `method == "L-BFGS-B"`
#' @return list with components:
#' - `par` final parameter vector
#' - `value` final log-likelihood
#' - `counts` number of function evaluations, gradient evaluations, and Hessian evaluations
#' - `convergence` convergence code
#' - `message` convergence message
#' - `hessian` estimated Hessian
#' @importFrom stats optim
#' @export
mle_lbfgsb_wei_series_md_c1_c2_c3 <- function(
    df,
    theta0,
    hessian = TRUE,
    ...,
    control = list()) {

    defaults <- list(fnscale = -1)
    control <- modifyList(defaults, control)

    optim(
        par = theta0,
        fn = function(theta) {
            loglik_wei_series_md_c1_c2_c3(df = df, theta = theta, ...)
        },
        gr = function(theta) {
            score_wei_series_md_c1_c2_c3(df = df, theta = theta, ...)
        },
        hessian = hessian,
        lower = rep(1e-3, length(theta0)),
        method = "L-BFGS-B",
        control = control)
}

#' Newton-raphson MLE solver for Weibull series model with masked component cause of
#' failure with observation of exact system lifetime and right censoring.
#' @param df data frame, right-censored lifetimes with masked component cause of failure
#' @param theta0 initial parameter vector
#' @param ... additional arguments passed to `loglik_wei_series_md_c1_c2_c3`,
#'            `score_wei_series_md_c1_c2_c3`, and `hessian_wei_series_md_c1_c2_c3`.
#' @param control list of control parameters to control the Newton solver.
#' @return list with components:
#' - `par` final parameter vector
#' - `value` final log-likelihood
#' - `iter` number of iterations
#' - `convergence` convergence code
#' - `hessian` estimated Hessian
#' @importFrom stats optim
#' @importFrom MASS ginv
#' @export
mle_newton_wei_series_md_c1_c2_c3 <- function(
    df,
    theta0,
    ...,
    control = list()) {

    defaults <- list(
        lr = 1,
        maxit = 100L,
        zero_tol = 1e-10,
        debug = FALSE,
        REPORT = 100L)
    control <- modifyList(defaults, control)

    convergence <- FALSE
    l <- loglik_wei_series_md_c1_c2_c3(df = df, theta = theta0, ...)
    iter <- 0L
    repeat {
        s <- score_wei_series_md_c1_c2_c3(df = df, theta = theta0, ...)
        h <- hessian_wei_series_md_c1_c2_c3(df = df, theta = theta0, ...)
        if (max(abs(s)) < control$zero_tol) {
            convergence <- TRUE
            break
        }
        if (iter > control$maxit) {
            break
        }    
        
        d <- ginv(h) %*% s
        a <- control$lr
        repeat {
            iter <- iter + 1L
            if (control$debug && iter %% control$REPORT == 0) {
                cat("iter:", iter, "l:", l, "theta:", theta0, ", s: ", s, "\n")
            }
            theta1 <- theta0 - a * d
            if (all(theta1 > 0)) {
                l1 <- loglik_wei_series_md_c1_c2_c3(df = df, theta = theta1, ...)
                if (l1 >= l) {
                    theta0 <- theta1
                    l <- l1
                    break
                }
            }
            if (iter > control$maxit) {
                break
            }
            a <- a / 2
        }
    }
    list(par = theta0, value = l, counts = iter,
         convergence = convergence, hessian = h, score = s)
}

#' Gradient ascent MLE solver for Weibull series model with masked component cause of
#' failure with observation of exact system lifetime and right censoring.
#' @param df data frame, right-censored lifetimes with masked component cause of failure
#' @param theta0 initial parameter vector
#' @param ... additional arguments passed to `loglik_wei_series_md_c1_c2_c3` and
#'            `score_wei_series_md_c1_c2_c3`.
#' @param control list of control parameters to control the gradient ascent solver.
#' @return list with components:
#' - `par` final parameter vector
#' - `value` final log-likelihood
#' - `iter` number of iterations
#' - `convergence` convergence code
#' - `hessian` estimated Hessian
#' @importFrom stats optim
#' @export
mle_grad_wei_series_md_c1_c2_c3 <- function(
    df,
    theta0,
    ...,
    control = list()) {

    defaults <- list(
        lr = 1,
        maxit = 100L,
        zero_tol = 1e-10,
        debug = FALSE,
        REPORT = 100L)

    control <- modifyList(defaults, control)
    convergence <- FALSE
    l <- loglik_wei_series_md_c1_c2_c3(df = df, theta = theta0, ...)
    iter <- 0L
    repeat {
        s <- score_wei_series_md_c1_c2_c3(df = df, theta = theta0, ...)
        if (max(abs(s)) < control$zero_tol) {
            convergence <- TRUE
            break
        }
        if (iter > control$maxit) {
            break
        }    
        a <- control$lr
        repeat {
            iter <- iter + 1L
            if (control$debug && iter %% control$REPORT == 0) {
                cat("iter:", iter, "l:", l, "theta:", theta0, ", s: ", s, "\n")
            }
            theta1 <- theta0 + a * s
            if (all(theta1 > 0)) {
                l1 <- loglik_wei_series_md_c1_c2_c3(df = df, theta = theta1, ...)
                if (l1 >= l) {
                    theta0 <- theta1
                    l <- l1
                    break
                }
            }
            a <- a / 2
        }
    }
    list(par = theta0, value = l, counts = iter, convergence = convergence,
        hessian = hessian_wei_series_md_c1_c2_c3(df = df, theta = theta0, ...),
        score = s)
}

#' Nelder-Mead MLE solver for the Weibull series model with
#' masked component cause of failure with observation of exact system lifetime and
#' right censoring.
#' @param df data frame, right-censored lifetimes with masked component cause of failure
#' @param theta0 initial parameter vector
#' @param ... additional arguments passed to `loglik_wei_series_md_c1_c2_c3`.
#' @param control list of control parameters passed to `optim` for `Nelder-Mead'
#' @return list with components:
#' - `par` final parameter vector
#' - `value` final log-likelihood
#' - `counts` number of function evaluations, gradient evaluations, and Hessian evaluations
#' - `convergence` convergence code
#' - `message` convergence message
#' - `hessian` estimated Hessian
#' @importFrom stats optim
#' @export
mle_nelder_wei_series_md_c1_c2_c3 <- function(
    df,
    theta0,
    ...,
    control = list()) {

    defaults <- list(fnscale = -1)
    control <- modifyList(defaults, control)

    optim(
        par = theta0,
        fn = function(theta) {
            loglik_wei_series_md_c1_c2_c3(df = df, theta = theta, ...)
        },
        hessian = TRUE,
        method = "Nelder-Mead",
        control = control)
}

#' Simulated annealing MLE solver for the Weibull series model with
#' masked component cause of failure with observation of exact system lifetime and
#' right censoring.
#' @param df data frame, right-censored lifetimes with masked component cause of failure
#' @param theta0 initial parameter vector
#' @param ... additional arguments passed to `loglik_wei_series_md_c1_c2_c3`.
#' @param control list of control parameters passed to `optim` for `SANN`
#' @return list with components:
#' - `par` final parameter vector
#' - `value` final log-likelihood
#' - `counts` number of function evaluations, gradient evaluations, and Hessian
#'   evaluations
#' - `convergence` convergence code
#' - `message` convergence message
#' - `hessian` estimated Hessian
#' @importFrom stats optim
#' @export
mle_sann_wei_series_md_c1_c2_c3 <- function(
    df,
    theta0,
    ...,
    control = list()) {

    defaults <- list(fnscale = -1)
    control <- modifyList(defaults, control)

    optim(
        par = theta0,
        fn = function(theta) {
            loglik_wei_series_md_c1_c2_c3(df = df, theta = theta, ...)
        },
        hessian = FALSE,
        method = "SANN",
        control = control)
}
