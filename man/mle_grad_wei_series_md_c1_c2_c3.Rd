% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/md_mle_wei_series_c1_c2_c3.R
\name{mle_grad_wei_series_md_c1_c2_c3}
\alias{mle_grad_wei_series_md_c1_c2_c3}
\title{Gradient ascent MLE solver for Weibull series model with masked component cause of
failure with observation of exact system lifetime and right censoring.}
\usage{
mle_grad_wei_series_md_c1_c2_c3(
  df,
  theta0,
  hessian = TRUE,
  ...,
  control = list()
)
}
\arguments{
\item{df}{data frame, right-censored lifetimes with masked component cause of failure}

\item{theta0}{initial parameter vector}

\item{hessian}{logical, if \code{TRUE}, then the Hessian is estimated and returned}

\item{...}{additional arguments passed to \code{loglik_wei_series_md_c1_c2_c3} and
\code{score_wei_series_md_c1_c2_c3}.}

\item{control}{list of control parameters to control the gradient ascent solver.}
}
\value{
list with components:
\itemize{
\item \code{par} final parameter vector
\item \code{value} final log-likelihood
\item \code{iter} number of iterations
\item \code{convergence} convergence code
\item \code{hessian} estimated Hessian
}
}
\description{
Gradient ascent MLE solver for Weibull series model with masked component cause of
failure with observation of exact system lifetime and right censoring.
}
