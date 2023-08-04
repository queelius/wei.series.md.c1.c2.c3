# This is the example series system used in the paper.
# It has five components in series. It is based on the
# example in Guo, table 2, but with two new components
# added. Each component has approximately the same
# reliability, so that there is no clear weakest link.
#
#   theta <- c(
#       k1 = 1.2576, lambda1 = 994.3661,
#       k2 = 1.1635, lambda2 = 908.9458,
#       k3 = 1.1308, lambda3 = 840.1141,
#       k4 = 1.1308, lambda4 = 940.1342,
#       k5 = 1.1308, lambda5 = 923.1631)

library(usethis)
alex_weibull_series <- list(
    shapes = c(1.2576, 1.1635, 1.1308, 1.1802, 1.2034),
    scales = c(994.3661, 908.9458, 840.1141, 940.1342, 923.1631),
    family = "weibull series system")

usethis::use_data(alex_weibull_series, overwrite = TRUE)
