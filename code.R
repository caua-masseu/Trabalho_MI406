library(robustbase)
library(MASS)

set.seed(10)
n <- 50
#Exemplo Simples
beta <- c(0.6, 2)
x1 <- rnorm(n = n, sd = 1, mean = 2)

size = 10

X <- matrix(c(rep(1, n), x1), nrow = n)
y <- X %*% beta + rnorm(n, sd = 1)

#y[n/5] <- size + y[n/5]
#y[floor(n/4)] <- size + y[n/4]
n_outliers <- 4
y[1:n_outliers] <- size + y[1:n_outliers]


MQO <- lm(y ~ X - 1)

#Detecção de Outlier e pontos influentes
which(abs(rstandard(MQO)) > 2)

plot(y ~ x1)
abline(a = coef(MQO)[1], b = coef(MQO)[2], col = "#FF0000")
L1 <-lmrob.lar(y = y, x = X)
abline(a = coef(L1)[1], b = coef(L1)[2], col = "#00FF00")
M <- rlm(y = y, x = X, psi = psi.huber)
abline(a = coef(M)[1], b = coef(M)[2], col = "#0000FF")
S <- lmrob.S(x = X, y = y, control = lmrob.control(psi = "bisquare"))
abline(a = coef(S)[1], b = coef(S)[2], col = "#9999FF")


size = 7.5

#Simulacao
n_rep <- 1000
bhatmqo <- data.frame()
bhatl1 <- data.frame()
bhatm <- data.frame()
bhats <- data.frame()
for(i in 1:n_rep){
  set.seed(i)
  # parâmetros
  beta <- c(1, 2.5, -1.5)
  # dados
  x1 <- rnorm(n = n, sd = 1, mean = 3)
  x2 <- rnorm(n = n, sd = 1, mean = 4)
  X <- matrix(c(rep(1, n), x1, x2), nrow = n)
  y <- X %*% beta + rnorm(n, sd = 1)
  # outliers
  y[(n-4):n] <- size + y[(n-4):n]
  # estimadores
  MQO <- lm(y ~ X - 1)
  L1 <-lmrob.lar(y = y, x = X)
  M <- rlm(y = y, x = X, psi = psi.huber)
  S <- lmrob.S(x = X, y = y, control = lmrob.control(psi = "bisquare"))
  bhatmqo <- rbind(bhatmqo, c(coef(MQO)[1:3]))
  bhatl1 <- rbind(bhatl1, c(coef(L1)[1:3]))
  bhatm <- rbind(bhatm, c(coef(M)[1:3]))
  bhats <- rbind(bhats, c(coef(S)[1:3]))
}
#colnames(bhatmqo) <- colnames(bhatl1) <- colnames(bhatm) <- colnames(bhats) <- c('beta_0$', '$\beta_0$', '$\beta_0$')
vies <- data.frame()
vies <- rbind(vies, apply(bhatmqo, MARGIN = 2, mean) - beta, 
              apply(bhatl1, MARGIN = 2, mean) - beta, 
              apply(bhatm, MARGIN = 2, mean) - beta, 
              apply(bhats, MARGIN = 2, mean) - beta)
colnames(vies) <- c('$beta_0$', '$beta_0$', '$beta_0$')
rownames(vies) <- c('MQO', 'L1', 'M', 'S')

xtable::xtable(vies)

plot(y ~x1)
hist(y)

lm(y ~ 1 + x1 + x2)
lmrob.lar(y = y, x = X)
MASS::rlm(y = y, x = X)
lmrob.S(x = X, y = y, control = lmrob.control("KS2011", max.it = 1000))[1]
lmrob(y ~ X - 1)



