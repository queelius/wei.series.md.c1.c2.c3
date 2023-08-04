# The source of this data set comes from Guo, table 2. When you use
# likelihood model that assumes a Weibull series system and candidate
# sets represented by Boolean vectors (x1, x2, x3) satisfy conditions
# C1, C2, and C3, the maximum likelihood estimates of the shape and
# scale parameters are:
#
#   theta.hat <- c(k1 = 1.2576, lambda1 = 994.3661,
#       k2 = 1.1635, lambda2 = 908.9458
#       k3 = 1.1308, lambda3 = 840.1141)
#
# This has a log-likelihood of -228.6851.

library(usethis)
library(readr)

guo_weibull_series_table_2 = read_csv("./data-raw/guo_weibull_series_table_2.csv")
usethis::use_data(guo_weibull_series_table_2, overwrite = TRUE)

guo_weibull_series_md <- list(
    mle = c(1.2576, 994.3661, 1.1635, 908.9458, 1.1308, 840.1141),
    shape_mles = c(1.2576, 1.1635, 1.1308),
    scale_mles = c(994.3661, 908.9458, 840.1141),
    loglike = -228.6851,
    p_hat = 0.215,
    tau = Inf,
    data = guo_weibull_series_table_2,
    family = "weibull series system")

usethis::use_data(guo_weibull_series_md, overwrite = TRUE)
