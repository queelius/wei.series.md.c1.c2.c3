#' Generates right-censored system failure times and right-censoring
#' indicators for a series system with the given data frame of
#' component lifetimes.
#'
#' @param df a data frame with the indicated component lifetimes
#' @param tau vector of right-censoring times, defaults to
#'            `NULL` (no right censoring)
#' @param control list of control parameters
#'  - `comp` - component lifetime prefix variable name, defaults to `t`, e.g.,
#'     `t1, t2, t3`.
#'  - `lifetime` - system lifetime variable name, defaults to `t`
#'  - `right_censoring_indicator` - right-censoring indicator variable, defaults
#'    to `delta`
#' @return (masked) data frame with masked data as described in the paper
#' @importFrom md.tools md_decode_matrix md_mark_latent
#' @importFrom dplyr %>%
#' @export
md_series_lifetime_right_censoring <- function(
    df,
    tau = NULL,
    control = list(comp = "t",
                   lifetime = "t",
                   right_censoring_indicator = "delta")) {
    # retrieve component lifetimes as a matrix
    t <- md_decode_matrix(df, control$comp)
    if (is.null(t)) {
        stop("comp not in colnames(df)")
    }

    s <- apply(t, 1, min)
    if (is.null(tau)) {
        df[, control$lifetime] <- s
    } else {
        df[, control$lifetime] <- ifelse(s < tau, s, tau)
        df[, control$right_censoring_indicator] <- s < tau
    }
    df %>% md_mark_latent(paste0(control$comp, 1:(ncol(t))))
}
