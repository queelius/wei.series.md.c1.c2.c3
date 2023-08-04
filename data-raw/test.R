theta.star <- c(shape1 = 1.2576, scale1 = 994.3661,
                shape2 = 1.1635, scale2 = 908.9458,
                shape3 = NA, scale3 = 500, #840.1141,
                shape4 = 1.1802, scale4 = 940.1342,
                shape5 = 1.2034, scale5 = 923.1631)

shapes <- theta.star[seq(1, length(theta.star), 2)]
scales <- theta.star[seq(2, length(theta.star), 2)]
df <- read.csv("./data-raw/data-shape-all.csv")
shapes3 <- unique(df$shape3)
tbl <- matrix(NA, nrow = length(shapes3), ncol = 4)
i <- 1L
for (shape3 in shapes3) {
    theta <- theta.star
    theta["shape3"] <- shape3
    shapes <- theta[seq(1, length(theta.star), 2)]
    scales <- theta[seq(2, length(theta.star), 2)]
    mttf.comp3 <- wei_mttf(shape = shapes[3], scale = scales[3])
    mttf.series <- wei_series_mttf(shapes = shapes, scales = scales)
    prob <- wei_series_cause(k = 3L, shapes = shapes, scales = scales)
    tbl[i, ] <- c(shape3, mttf, mttf.series, prob)
    i <- i + 1L
}
colnames(tbl) <- c("shape[3]", "mttf[3]", "mttf", "prob[3]")
print(tbl)






















### for scale vary

theta.star <- c(shape1 = 1.2576, scale1 = 994.3661,
                shape2 = 1.1635, scale2 = 908.9458,
                shape3 = 1.1308, scale3 = NA,
                shape4 = 1.1802, scale4 = 940.1342,
                shape5 = 1.2034, scale5 = 923.1631)

shapes <- theta.star[seq(1, length(theta.star), 2)]
scales <- theta.star[seq(2, length(theta.star), 2)]
df <- read.csv("./data-raw/data-shape-all.csv")
shapes3 <- unique(df$shape3)
tbl <- matrix(NA, nrow = length(shapes3), ncol = 4)
i <- 1L
for (shape3 in shapes3) {
    theta <- theta.star
    theta["shape3"] <- shape3
    shapes <- theta[seq(1, length(theta.star), 2)]
    scales <- theta[seq(2, length(theta.star), 2)]
    mttf.comp3 <- wei_mttf(shape = shapes[3], scale = scales[3])
    mttf.series <- wei_series_mttf(shapes = shapes, scales = scales)
    prob <- wei_series_cause(k = 3L, shapes = shapes, scales = scales)
    tbl[i, ] <- c(shape3, mttf, mttf.series, prob)
    i <- i + 1L
}
colnames(tbl) <- c("shape[3]", "mttf[3]", "mttf", "prob[3]")
print(tbl)

