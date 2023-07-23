#' Bernoulli candidate set model is a particular type of *uninformed* model.
#' Note that we do not generate candidate sets with this function. See
#' `md_cand_sampler` for that.
#'
#' This model satisfies conditions C1, C2, and C3.
#' The failed component will be in the corresponding candidate set with
#' probability 1, and the remaining components will be in the candidate set
#' with probability `p`.
#'
#' @param df masked data.
#' @param p a vector of probabilities (p[j] is the probability that the jth
#'         system will include a non-failed component in its candidate set,
#'         assuming the jth system is not right-censored).
#' @param control 
#'  - `prob` - column prefix for component probabilities, defaults to
#'    `q`, e.g., `q1, q2, q3`.
#'  - `comp` - column prefix for component lifetimes, defaults to `t`,
#'     e.g., `t1, t2, t3`.
#'  - `right_censoring_indicator` right-censoring indicator column name, if
#'     TRUE, then the system lifetime is right-censored, otherwise it is
#'     observed.
#' @importFrom md.tools md_decode_matrix md_encode_matrix md_mark_latent
#' @importFrom dplyr %>% bind_cols
#' @export
md_bernoulli_cand_c1_c2_c3 <- function(
    df,
    p,
    control = list()) {

    defaults <- list(
        comp = "t",
        prob = "q",
        right_censoring_indicator = NULL)
    control <- modifyList(defaults, control)

    n <- nrow(df)
    if (n == 0) {
        return(df)
    }
    p <- rep(p, length.out = n)
    Tm <- md_decode_matrix(df, control$comp)
    if (is.null(Tm)) {
        stop("No component lifetime variables found")
    }
    m <- ncol(Tm)
    Q <- matrix(p, nrow = n, ncol = m)
    Q[cbind(1:n, apply(Tm, 1, which.min))] <- 1

    if (!is.null(control$right_censoring_indicator)) {

        if (!control$right_censoring_indicator %in% colnames(df)) {
            stop("right_censoring_indicator variable not in colnames(df)")
        }
        delta <- df[[control$right_censoring_indicator]]
        Q[delta, ] <- 0
    }

    # remove in case it already has columns for q1,...,qm
    df[ , paste0(control$prob, 1:m)] <- NULL
    df %>% bind_cols(md_encode_matrix(Q, control$prob)) %>%
           md_mark_latent(paste0(control$prob, 1:m))
}

#' Candidate set generator. Requires columns for component probabilities
#' e.g., `q1,...,qm` where `qj` is the probability that the jth component
#' will be in the corresponding candidate set generated for that observation
#' in the `md` table.
#'
#' @param df (masked) data frame
#' @param control 
#'  - `prob` - column prefix for component probabilities, defaults to
#'    `q`, e.g., `q1, q2, q3`.
#'  - `candset` - column prefix for candidate sets (as Boolean matrix),
#'    defaults to `x`, e.g., `x1, x2, x3`.
#' @importFrom md.tools md_decode_matrix md_encode_matrix
#' @importFrom dplyr %>% bind_cols
#' @importFrom stats runif
#' @export
md_cand_sampler <- function(df, control = list()) {

    if (!is.data.frame(df)) {
        stop("df must be a data frame")
    }

    n <- nrow(df)
    if (n == 0) {
        return(df)
    }

    defaults <- list(
        prob = "q",
        candset = "x")
    control <- modifyList(defaults, control)

    Q <- md_decode_matrix(df, control$prob)
    m <- ncol(Q)
    if (m == 0) {
        stop("No component probabilities found")
    }

    X <- matrix(NA, nrow=n, ncol=m)
    for (i in 1:n)
        X[i, ] <- runif(m) <= Q[i, ]
    df %>% bind_cols(md_encode_matrix(X, control$candset))
}
