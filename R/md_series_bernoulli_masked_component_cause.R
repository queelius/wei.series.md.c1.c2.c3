#' Bernoulli candidate set model is a particular type of *uninformed* model.
#' Note that we do not generate candidate sets with this function. See
#' `md_cand_sampler` for that.
#'
#' This model satisfies conditions C1, C2, and C3.
#' The failed component will be in the corresponding candidate set with
#' probability 1, and the remaining components will be in the candidate set
#' with probability `p` (the same probability for each component). `p` 
#' may be different for each system, but it is assumed to be the same for
#' each component within a system, so `p` can be a vector such that the
#' length of `p` is the number of systems in the data set (with recycling
#' if necessary).
#'
#' @param df masked data.
#' @param p a vector of probabilities (p[j] is the probability that the jth
#'          system will include a non-failed component in its candidate set,
#'          assuming the jth system is not right-censored).
#' @param prob column prefix for component probabilities, defaults to
#'             `q`, e.g., `q1, q2, q3`.
#' @param comp column prefix for component lifetimes, defaults to `t`,
#'             e.g., `t1, t2, t3`.
#' @param right_censoring_indicator right-censoring indicator column name.
#'     if TRUE, then the system lifetime is right-censored, otherwise it is
#'     observed. If NULL, then no right-censoring is assumed. Defaults to
#'     `delta`.
#' @importFrom md.tools md_decode_matrix md_encode_matrix md_mark_latent
#' @importFrom dplyr %>% bind_cols
#' @export
md_bernoulli_cand_c1_c2_c3 <- function(
    df,
    p,
    prob = "q",
    comp = "t",
    right_censoring_indicator = "delta") {

    n <- nrow(df)
    if (n == 0) {
        return(df)
    }
    p <- rep(p, length.out = n)
    Tm <- md_decode_matrix(df, comp)
    if (is.null(Tm)) {
        stop("No component lifetime variables found")
    }
    m <- ncol(Tm)
    Q <- matrix(p, nrow = n, ncol = m)
    Q[cbind(1:n, apply(Tm, 1, which.min))] <- 1

    if (!is.null(right_censoring_indicator)) {

        if (!right_censoring_indicator %in% colnames(df)) {
            stop("right_censoring_indicator not in colnames(df)")
        }
        Q[!df[[right_censoring_indicator]], ] <- 0
    }

    # remove in case it already has columns for q1,...,qm
    df[ , paste0(prob, 1:m)] <- NULL
    df %>% bind_cols(md_encode_matrix(Q, prob)) %>%
           md_mark_latent(paste0(prob, 1:m))
}

#' Candidate set generator. Requires columns for component probabilities
#' e.g., `q1,...,qm` where `qj` is the probability that the jth component
#' will be in the corresponding candidate set generated for that observation
#' in the `md` table.
#'
#' @param df (masked) data frame
#' @param prob column prefix for component probabilities, defaults to
#'             `q`, e.g., `q1, q2, q3`.
#' @param candset column prefix for candidate sets (as Boolean matrix),
#'                defaults to `x`, e.g., `x1, x2, x3`.
#' @importFrom md.tools md_decode_matrix md_encode_matrix
#' @importFrom dplyr %>% bind_cols
#' @importFrom stats runif
#' @export
md_cand_sampler <- function(df, prob = "q", candset = "x") {

    if (!is.data.frame(df)) {
        stop("df must be a data frame")
    }

    n <- nrow(df)
    if (n == 0) {
        return(df)
    }

    Q <- md_decode_matrix(df, prob)
    m <- ncol(Q)
    if (m == 0) {
        stop("No component probabilities found")
    }

    X <- matrix(NA, nrow=n, ncol=m)
    for (i in 1:n) {
        X[i, ] <- runif(m) <= Q[i, ]
    }
    df %>% bind_cols(md_encode_matrix(X, candset))
}

#' Check if a masked data frame is identifiable.
#' 
#' We want to make sure that the right-censored, masked component cause of failure
#' data has enough information in it for the likelihood model to be uniquely
#' identified.
#'
#' @param df masked data frame
#' @param candset column prefix for candidate sets, defaults to `x`,
#' e.g., `x1, x2, x3`.
#' @return returns a string describing the data as yielding an identifiable
#'         likelihood function, or one or more reasons why it is not identifiable.
#' @importFrom md.tools md_decode_matrix
#' @export
md_bernoulli_cand_c1_c2_c3_identifiable <- function(df, candset = "x") {

    C <- md_decode_matrix(df, candset)
    if (is.null(C)) {
        stop("No candidate set variables found")
    }
    n <- nrow(C)
    if (n == 0) {
        return("No observations")
    }

    reasons <- c()
    # make sure a component appears at least once in
    # a candidat set by checking each column of C
    # and verifying that it has at least one 1 (TRUE)
    # this also checks to make sure that not every
    # observation is right-censored, since if that
    # were the case, then every column of C would
    # be all 0 (FALSE)
    if (any(colSums(C) == 0)) {
        reasons <- c(reasons, "No component appears in a candidate set")
    }

    # make sure that there is at least one observation
    # (row) that has a 0 (FALSE) in a column of C
    # if every candidate set has all components,
    # then the candidate sets convey no information
    # about the components, and the model is not
    # identifiable.
    if (all(rowSums(C) == ncol(C))) {
        reasons <- c(reasons, "Every candidate set has all components")
    }


    # if a column in C is all 1 (TRUE), then the
    # corresponding component is always in the
    # candidate set, and the model is not
    # identifiable.
    cs <- colSums(C) == nrow(C)
    if (any(cs)) {
        reasons <- c(reasons, paste0("Component(s) ", paste(which(cs), collapse = ", "), " always in candidate set"))
    }

    # these conditions are necessary, but not
    # sufficient for identifiability. it's not clear
    # that we can determine ahead of time when
    # identifiability is not possible, but if the
    # MLE is not unique, that also means that the
    # likelihood function is not identifiable.
    # we do not check for that here, since it
    # depends on the likelihood model, e.g., weibull series,
    # and so on. at the time we find the MLE, though,
    # we can check to see if the likelihood function
    # is approximately flat where it is at a maximum,
    # and if so, then the likelihood function is not identifiable.
    return(reasons)
}


#' Add info about component cause of failure to a
#' masked data frame with latent components available.
#' @param df masked data frame
#' @param comp column prefix for component lifetimes, defaults to `t`,
#' @return a data frame with a new `k` column
#' @importFrom md.tools md_decode_matrix
#' @importFrom dplyr %>% bind_cols
#' @export
md_decorate_component_cause_k <- function(df, comp = "t") {

    Tm <- md_decode_matrix(df, comp)
    if (is.null(Tm)) {
        stop("No component lifetime variables found")
    }
    K <- apply(Tm, 1, which.min)
    df %>% bind_cols(k = K) %>% md_mark_latent("k")
}

#' Empirical sampler for C[i] | K[i] = k given a masked
#' data frame with component cause of failures shown.
#' 
#' In the likelihood model, we assume rather that
#' C[i] | T[i] = t, K[i] = k satisfies certain assumptions.
#' For tractability and efficiency, we're discarding the conditional relation on
#' T[i] and only keeping the conditional relation on K[i] by using the
#' empirical distribution of C[i] | K[i] = k.
#' 
#' If k is not known, then we can sample from the empirical distribution of
#' C[i], which may still be reasonable, since by Condition 2,
#' Pr{C[i] = c[i] | T[i] = t[i], K[i] = j} =
#' Pr{C[i] = c[i] | T[i] = t[i], K[i] = j'}  for all j, j' in c[i],
#' and so in general varying the component cause among the components
#' in the candidate set does not change the probability of the candidate set.
#' 
#' If the model is more informed, then dropping the conditional relation on
#' K[i] will lose information, and may bias the results in complex ways.
#' 
#' @param df masked data frame
#' @param candset column prefix for candidate sets, defaults to `x`,
#'                e.g., `x1, x2, x3`.
#' @param cause column name for component cause of failure, defaults to `k`
#' @return a data frame of candidate sets sampled from the empirical
#'         distribution of C[i] | K[i] = k
#' @importFrom md.tools md_decode_matrix
#' @importFrom dplyr %>% starts_with select sample_n
#' @export
md_sample_candidates <- function(df, n, k = NA, cause = "k", candset = "x") {
    df[df[[cause]] == k,] %>% select(starts_with(candset)) %>%
        sample_n(size = n, replace = TRUE)
}

#' Conditional sampler C[i] | T[i] = t, K[i] = k
#' We have a smoothing parameter that specifies the bin's width for the
#' discretization of the system lifetimes.
#' 
#' If component cause of failure is known, then we can sample from the
#' conditional distribution of C[i] | T[i] = t, K[i] = k, otherwise
#' we can sample from the empirical distribution of C[i] | T[i] = t[i],
#' which may still be reasonable, since by Condition 2,
#' Pr{C[i] = c[i] | T[i] = t[i], K[i] = j} = 
#' Pr{C[i] = c[i] | T[i] = t[i], K[i] = j'}  for all j, j' in c[i].
#' Of course, we may violate Condition 1, Pr{K[i] in C[i]} = 1.
#' 
#' Note that if |c[i]| = 1, then we can sample from the empirical
#' distribution of C[i] | T[i] = t[i], K[i] = j, since there is only
#' one possible value for K[i], K[i] = j if c[i] = {j}. Otherwise,
#' if |c[i]| > 1, and we don't know the component cause of failure,
#' then we can sample from the empirical distribution of C[i] | T[i] = t[i].
#' 
#' In the semi-parametric bootstrap, we generate the samples ourselves
#' from our estimate of theta, so we can sample from the conditional
#' distribution of C[i] | T[i] = t[i], K[i] = j, since we know the
#' simulated component cause of failure. This is the main reason we
#' want to sample C[i] | T[i] = t[i], K[i] = j, so this is fine.
#'
#' @param t observed lifetime, defaults to NA (unknown or unconditional)
#' @param df data frame (sample) that we used to estimate C[i] | T[i] = t[i], K[i] = j
#'           `df` should only contain data in which the system failed, rather than being
#'           right-censored.
#' @param k component cause of failure, defaults to NA (unknown)
#' @param nbins number of bins to use for discretizing the component lifetimes
conditional_masked_cause <- function(n, df, t, k, nbins = 10) {

  df <- df[df$delta, ] %>% mutate(
    bins = cut(t, breaks = quantile(t, probs = seq(0, 1, 1/nbins)),
        include.lowest = TRUE))
  
  # Find which bin t belongs to
  bin <- df$bins[df$t == t]
  
  # Sample indices to so that we know which rows in `df` to sample from
  indices <- sample(df$bins == bin & df$k == k, n)

  # Sample from the empirical distribution of C[i] | T[i] = t[i], K[i] = j
    df[indices, ] %>% select(starts_with("x")) %>% sample_n(n, replace = TRUE)
}
