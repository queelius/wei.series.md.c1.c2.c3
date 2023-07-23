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
#'   \item{t}{Time of failure of series system}
#'   \item{x1}{Indicates whether component 1 is in the candidate set}
#'   \item{x2}{Indicates whether component 2 is in the candidate set}
#'   \item{x3}{Indicates whether component 3 is in the candidate set}
#' }
#' 
#' @seealso \code{\link{guo_weibull_series_mle}}
#'
#' @source H. Guo, F. Szidarovszky, and P. Niu,
#' "Estimating component reliabilities from incomplete system failure data,"
#' in 2013 Proceedings Annual Reliability and Maintainability Symposium (RAMS),
#' 2013. [Online]. Available: https://doi.org/10.1109/rams.2013.6517765
#'
#' @examples
#' head(guo_weibull_series_md)
#' loglik <- md_loglike_weibull_series_C1_C2_C3(guo_weibull_series_md
#'     deltavar = NULL) # no right-censoring in this data set
#' sol <- optim(par = guo_weibull_series_mle$mle,
#'              fn = loglik,
#'              hessian = TRUE,
#'              control = list(parscale = c(1, 1000, 1, 1000, 1, 1000)))
#' abs(sol$value - guo_weibull_series_mle$loglike) < 1e-4
#' abs(sol$par - guo_weibull_series_mle$mle) < 1e-4
"guo_weibull_series_md"

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
#'  \item{family}{The family of the likelihood model}
#' }
#' @seealso \code{\link{guo_weibull_series_md}}
"guo_weibull_series_mle"