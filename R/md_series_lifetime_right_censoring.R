#' Generates right-censored system failure times and right-censoring
#' indicators for a series system with the given data frame of
#' component lifetimes.
#'
#' @param df a data frame with the indicated component lifetimes
#' @param tau vector of right-censoring times, defaults to
#'            `Inf` (no right censoring)
#' @param comp component lifetime prefix variable name, defaults to `t`, e.g.,
#'             `t1, t2, t3`.
#' @param lifetime system lifetime variable name, defaults to `t`
#' @param right_censoring_indicator right-censoring indicator variable, defaults
#'                                  to `delta`
#' @return (masked) data frame with masked data as described in the paper
#' @importFrom md.tools md_decode_matrix md_mark_latent
#' @importFrom dplyr %>%
#' @export
md_series_lifetime_right_censoring <- function(
    df,
    tau = Inf,
    comp = "t",
    lifetime = "t",
    right_censoring_indicator = "delta") {

    t <- md_decode_matrix(df, comp)
    if (is.null(t)) {
        stop("comp not in colnames(df)")
    }

    s <- apply(t, 1, min)
    df[, lifetime] <- ifelse(s < tau, s, tau)
    df[, right_censoring_indicator] <- s < tau
    df %>%
        md_mark_latent(paste0(comp, 1:(ncol(t)))) %>%
        md_mark_latent(lifetime)
}
