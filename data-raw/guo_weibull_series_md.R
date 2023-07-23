# The source of this data set comes from Guo, table 2. When you use
# likelihood model that assumes a Weibull series system and candidate
# sets represented by Boolean vectors (X1, X2, X3) satisfy conditions
# C1, C2, and C3, the maximum likelihood estimates of the shape and
# scale parameters are:
#   - β1 = 1.2576 η1 = 994.3661
#   - β2 = 1.1635 η2 = 908.9458
#   - β3 = 1.1308 η3 = 840.1141
#
#   - guo.theta <- c(1.2576, 994.3661, 1.1635, 908.9458, 1.1308, 840.1141)
#
# This has a log-likelihood of -228.6851.

library(usethis)
library(readr)

guo_weibull_series_md <- read_csv("./inst/guo_weibull_series_md.csv")
usethis::use_data(guo_weibull_series_md, overwrite = TRUE)


guo_weibull_series_mle <- list(
    mle = c(1.2576, 994.3661, 1.1635, 908.9458, 1.1308, 840.1141),
    shape_mles = c(1.2576, 1.1635, 1.1308),
    scale_mles = c(994.3661, 908.9458, 840.1141),
    loglike = -228.6851,
    p.hat = 0.215,
    data = guo_weibull_series_md,
    family = "weibull series system")
usethis::use_data(guo_weibull_series_mle, overwrite = TRUE)
