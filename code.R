library(robustbase)

set.seed(1)
n <- 100

beta <- c(0.6, 2, 3)
x1 <- rnorm(n = n, sd = 1)
x2 <- rnorm(n = n, sd = 1)

#Put a outlier in x_1 and x_2 in different times
size = 10
#x1[n/2] <- x1[n/2] + size
#x2[3*n/4] <- x2[3*n/4] + size

X <- matrix(c(rep(1, n), x1, x2), nrow = n)
y <- X %*% beta + rnorm(n)

y[n/5] <- 20*y[n/5]
y[n/4] <- 20*y[n/4]


plot(y ~x1)
hist(y)

lm(y ~ 1 + x1 + x2)
coef(lmrob.lar(y = y, x = X))
MASS::rlm(y = y, x = X)
lmrob.S(x = X, y = y, control = lmrob.control("KS2011", max.it = 1000))[1]
