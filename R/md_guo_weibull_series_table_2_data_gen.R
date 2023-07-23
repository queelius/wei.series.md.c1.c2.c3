#' Generate data like that from the Guo et al. model, table 2.
#' 
#' @param n Number of observations to generate. Default is 30, same as the
#'          Guo et al. data set.
#' @param shapes The shape parameter values. Default is the MLE from the Guo et al.
#'             data set, `1.2576, 1.1635, 1.1308`
#' @param scales The scale parameter values. Default is the MLE from the Guo et
#'               al. data set, `994.3661, 908.9458, 840.1141`
#' @param p The probability of a non-failed component being in the candidate set.
#'         Default is the estimated probability from the Guo et al. data set, `0.215`.
#' @param tau The right-censoring time. Default is `Inf`, same as the Guo et al.
#'            data set, which had no right-censoring.
#' @return A data frame
#' @importFrom stats rweibull
#' @importFrom md.tools md_encode_matrix
#' @export
generate_guo_weibull_table_2_data <- function(
    n = 30,
    shapes = guo_weibull_series_mle$shape_mles,
    scales = guo_weibull_series_mle$scale_mles,
    p = guo_weibull_series_mle$p.hat,
    tau = Inf) {

    m <- length(shapes)
    comp_times <- matrix(nrow = n, ncol = m)

    for (j in 1:m)
        comp_times[, j] <- rweibull(
            n = n,
            shape = shapes[j],
            scale = scales[j])
    comp_times <- md_encode_matrix(comp_times, "t")

    comp_times %>% md_series_lifetime_right_censoring(tau) %>%
        md_bernoulli_cand_c1_c2_c3(p) %>% md_cand_sampler()
}