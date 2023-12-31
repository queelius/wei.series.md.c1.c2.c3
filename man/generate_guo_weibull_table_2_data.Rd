% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/md_guo_weibull_series_table_2_data_gen.R
\name{generate_guo_weibull_table_2_data}
\alias{generate_guo_weibull_table_2_data}
\title{Generate data like that from the Guo et al. model, table 2.}
\usage{
generate_guo_weibull_table_2_data(
  n = nrow(guo_weibull_series_md$data),
  shapes = guo_weibull_series_md$shape_mles,
  scales = guo_weibull_series_md$scale_mles,
  p = guo_weibull_series_md$p_hat,
  tau = guo_weibull_series_md$tau
)
}
\arguments{
\item{n}{Number of observations to generate. Default is 30, same as the
Guo et al. data set in table 2.}

\item{shapes}{The shape parameter values. Default is the MLE from the Guo et al.
data set in table 2, \verb{1.2576, 1.1635, 1.1308}}

\item{scales}{The scale parameter values. Default is the MLE from the Guo et
al. data set in table 2, \verb{994.3661, 908.9458, 840.1141}}

\item{p}{The probability of a non-failed component being in the candidate set.
Default is the estimated probability from the Guo et al. data set
in Table 2, \code{0.215}.}

\item{tau}{The right-censoring time. Default is \code{Inf}, same as the Guo et al.
data set in Table 2, which had no right-censoring.}
}
\value{
A data frame
}
\description{
Generate data like that from the Guo et al. model, table 2.
}
