% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/md_series_bernoulli_masked_component_cause.R
\name{conditional_masked_cause}
\alias{conditional_masked_cause}
\title{Conditional candidate set sampler when component cause of failure is known
and the system lifetime is observed. (UNTESTED)}
\usage{
conditional_masked_cause(n, df, t, k, nbins = 10)
}
\arguments{
\item{df}{data frame (sample) that we used to estimate C\link{i} | T\link{i} = t\link{i}, K\link{i} = j
\code{df} should only contain data in which the system failed, rather than being
right-censored.}

\item{t}{observed lifetime, defaults to NA (unknown or unconditional)}

\item{k}{component cause of failure, defaults to NA (unknown)}

\item{nbins}{number of bins to use for discretizing the component lifetimes}
}
\description{
Conditional sampler C\link{i} | T\link{i} = t, K\link{i} = k
We have a smoothing parameter that specifies the bin's width for the
discretization of the system lifetimes.
}
\details{
If component cause of failure is known, then we can sample from the
conditional distribution of C\link{i} | T\link{i} = t, K\link{i} = k, otherwise
we can sample from the empirical distribution of C\link{i} | T\link{i} = t\link{i},
which may still be reasonable, since by Condition 2,
Pr{C\link{i} = c\link{i} | T\link{i} = t\link{i}, K\link{i} = j} =
Pr{C\link{i} = c\link{i} | T\link{i} = t\link{i}, K\link{i} = j'}  for all j, j' in c\link{i}.
Of course, we may violate Condition 1, Pr{K\link{i} in C\link{i}} = 1.

Note that if |c\link{i}| = 1, then we can sample from the empirical
distribution of C\link{i} | T\link{i} = t\link{i}, K\link{i} = j, since there is only
one possible value for K\link{i}, K\link{i} = j if c\link{i} = {j}. Otherwise,
if |c\link{i}| > 1, and we don't know the component cause of failure,
then we can sample from the empirical distribution of C\link{i} | T\link{i} = t\link{i}.

In the semi-parametric bootstrap, we generate the samples ourselves
from our estimate of theta, so we can sample from the conditional
distribution of C\link{i} | T\link{i} = t\link{i}, K\link{i} = j, since we know the
simulated component cause of failure. This is the main reason we
want to sample C\link{i} | T\link{i} = t\link{i}, K\link{i} = j, so this is fine.
}
