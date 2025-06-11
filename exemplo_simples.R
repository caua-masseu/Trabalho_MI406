library(robustbase)
library(MASS)

#Exemplo Simples
set.seed(10)
n <- 50
p <- 2
beta <- c(0.6, 2)
x1 <- rnorm(n = n, sd = 1, mean = 2)

size = 10

X <- matrix(c(rep(1, n), x1), nrow = n)
y <- X %*% beta + rnorm(n, sd = 1) #Desvio padrado de y = 1                     

n_outliers <- 4
y[11] <- size + y[11]
y[12] <- size + y[12]
y[41] <- size + y[41]
y[50] <- size + y[50]

MQO <- lm(y ~ X - 1)

which()

#Detecção de Outlier e pontos influentes

data_res <- data.frame(res =  residuals(MQO), r_std = rstandard(MQO), r_std_d = rstudent(MQO))[c(10,11,12,40,41,42,49,50),]
xtable::xtable(data_res, digits = 3)
which(abs(rstandard(MQO)) > 2)
which(abs(rstudent(MQO)) >  2)
cooks.distance(MQO)

point_inf <- data.frame(dffits = dffits(MQO), dfbeta0 = dfbetas(MQO)[,1], 
                        dfbeta1 = dfbetas(MQO)[,2], dist_cook = cooks.distance(MQO))[c(10,11,12,40,41,42,49,50),]

xtable::xtable(point_inf, digits = 3)

which(cooks.distance(MQO) > qf(p = 0.5, df1 = p, df2 = n - p))
which(dffits(MQO) > 2*sqrt(p/n))
which(dfbetas(MQO)[,1] > 2/sqrt(n)) #beta0
which(dfbetas(MQO)[,2] > 2/sqrt(n)) #beta1
inf <- influence.measures(MQO, )
summ <- summary(inf)

png(filename = "plot_exemplo_simples.png", width = 900, height = 600, )

plot(y ~ x1)
abline(a = coef(MQO)[1], b = coef(MQO)[2], col = "#FF0000")
L1 <-lmrob.lar(y = y, x = X)
abline(a = coef(L1)[1], b = coef(L1)[2], col = "#00FF00")
M <- rlm(y = y, x = X, psi = psi.huber)
abline(a = coef(M)[1], b = coef(M)[2], col = "#0000FF")
S <- lmrob.S(x = X, y = y, control = lmrob.control(psi = "bisquare"))
abline(a = coef(S)[1], b = coef(S)[2], col = "#9900FF")
legend("topleft", legend = c("MQO", "L1", "M", "S"),  fill = c('#FF0000', "#00FF00", "#0000FF", "#9900FF"), cex = 1)

dev.off()

xtable::xtable(data.frame(MQO = coef(MQO), L1 = coef(L1), M = coef(M), S = coef(S), row.names = c("$beta0$", "$beta1$")), digits = 3)
