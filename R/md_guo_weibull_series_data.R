#' Weibull series MLE of `guo_weibull_series_md`
#' 
#' When you use likelihood model that assumes a Weibull series system and 
#' candidate sets represented by Boolean vectors (x1, x2, x3) that satisfy 
#' conditions C1, C2, and C3, the MLE of the shape and
#' scale parameters are:
#'
#'   - β1 = 1.2576
#'   - η1 = 994.3661
#'   - β2 = 1.1635
#'   - η2 = 908.9458
#'   - β3 = 1.1308
#'   - η3 = 840.1141
#'
#'   - θ̂  = (β1, η1, β2, η2, β3, η3)
#'
#' This has a log-likelihood of -228.6851.
#' 
#' @format A list with the following components:
#' \describe{
#'  \item{mle}{The MLE of the shape and scale parameters}
#'  \item{loglike}{The log-likelihood of the MLE}
#'  \item{data}{The data used to compute the MLE}
#'  \item{p_hat}{The estimated probability of a non-failed component being in the candidate set}
#'  \item{tau}{The right-censoring time (infinity)}
#'  \item{family}{The family of the likelihood model (weibull series)}
#' }
#' @seealso \code{\link{guo_weibull_series_table_2}}
"guo_weibull_series_md"

#' Example Data for a Series System (Table 2, Guo 2013)
#' 
#' The source of this data set comes from "Estimating Component Reliabilities
#' from Incomplete System Failure Data", Table 2: Example Data for a Series
#' System.
#'
#' When you use likelihood model that assumes a Weibull series system and 
#' candidate sets represented by Boolean vectors (x1, x2, x3) that satisfy 
#' conditions C1, C2, and C3, the MLE of the shape and
#' scale parameters are:
#'
#'   - β1 = 1.2576
#'   - η1 = 994.3661
#'   - β2 = 1.1635
#'   - η2 = 908.9458
#'   - β3 = 1.1308
#'   - η3 = 840.1141
#'
#'   - θ̂  = (β1, η1, β2, η2, β3, η3)
#'
#' This has a log-likelihood of -228.6851.
#'
#' @format A data frame with 30 rows and 4 columns:
#' \describe{
#'   \item{t}{number, lifetime of series system}
#'   \item{x1}{logical, TRUE indicates component 1 is in the candidate set, FALSE otherwise}
#'   \item{x2}{logical, TRUE indicates component 2 is in the candidate set, FALSE otherwise}
#'   \item{x3}{logical, TRUE indicates component 3 is in the candidate set, FALSE otherwise}
#'   \item{delta}{logical, Right-censoring indicator, TRUE if the system is observed}
#' }
#' 
#' @seealso \code{\link{guo_weibull_series_md}}
#'
#' @source H. Guo, F. Szidarovszky, and P. Niu,
#' "Estimating component reliabilities from incomplete system failure data,"
#' in 2013 Proceedings Annual Reliability and Maintainability Symposium (RAMS),
#' 2013. [Online]. Available: https://doi.org/10.1109/rams.2013.6517765
#'
#' @examples
#' head(guo_weibull_series_md)
#' sol <- optim(par = guo_weibull_series_md$mle,
#'              fn = loglik_wei_series_md_c1_c2_c3,
#'              hessian = TRUE,
#'              control = list(fnscale = -1,
#'                             parscale = c(1, 1000, 1, 1000, 1, 1000)),
#'              df = guo_weibull_series_md$data)
#' abs(sol$value - guo_weibull_series_md$loglike) < 1e-4
#' abs(sol$par - guo_weibull_series_md$mle) < 1e-4
"guo_weibull_series_table_2"

#' This is the example series system used in the paper.
#' It has five components in series. It is based on the
#' example in Guo, table 2, but with two new components
#' added. Each component has approximately the same
#' reliability, so that there is no clear weakest link.
#' 
#' theta = (shape1 = 1.2576, scale1 = 994.3661,
#'          shape2 = 1.1635, scale2 = 908.9458,
#'          shape3 = 1.1308, scale3 = 840.1141,
#'          shape4 = 1.1802, scale4 = 940.1342,
#'          shape5 = 1.2034, scale5 = 923.1631)
#' 
#' @format A list with the following components:
#' \describe{
#'   \item{theta}{The shape and scale parameters of the components combined}
#'   \item{shapes}{The shape parameters of the components}
#'   \item{scales}{The scale parameters of the components}
#'   \item{family}{The family of the likelihood model (weibull series)}
#' }
#' @seealso \code{\link{guo_weibull_series_md}}
"alex_weibull_series"
