# This is the example series system used in the paper.
# It has five components in series. It is based on the
# example in Guo, table 2, but with two new components
# added. Each component has approximately the same
# reliability, so that there is no clear weakest link.
library(usethis)

theta <- c(shape1 = 1.2576, scale1 = 994.3661,
           shape2 = 1.1635, scale2 = 908.9458,
           shape3 = 1.1308, scale3 = 840.1141,
           shape4 = 1.1802, scale4 = 940.1342,
           shape5 = 1.2034, scale5 = 923.1631)
shapes <- theta[seq(1, length(theta), 2)]
scales <- theta[seq(2, length(theta), 2)]

alex_weibull_series <- list(
    theta = theta,
    shapes = shapes,
    scales = scales,
    family = "weibull series system")

usethis::use_data(alex_weibull_series, overwrite = TRUE)
