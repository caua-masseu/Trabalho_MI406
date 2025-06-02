set.seed(1)

n <- 300

x1 <- rnorm(n)
x2 <- rnorm(n)
beta <- c(0.5, 1, 2)

y <- beta[1] + beta[2]*x1 + beta[3]*x2 + rnorm(n)

#Put outlier in the model
size <- 7.5

x1[n/2] <- x1[n/2] + size
x2[2*n/3] <- x2[2*n/3] + size

y <- beta[1] + beta[2]*x1 + beta[3]*x2 + rnorm(n)

lm(y ~ 1 + x1 + x2)
