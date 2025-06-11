data("LifeCycleSavings")

library(ggplot2)
library(tidyverse)
library(GGally)

library(robustbase)
library(MASS)

# Analise exploratoria

corrplot::corrplot(cor(LifeCycleSavings[,]))
ggpairs(LifeCycleSavings[, -c(3)]) 
lcs <- LifeCycleSavings[, -c(3)]

# MQO 
ajuste <- lm(sr ~ pop15 + dpi + ddpi, data = lcs)
paises <- row.names(lcs)

# Detecção de outliers

outlier <- function(fit, t_critical = 2, alpha = 0.05){
  n <- length(residuals(fit)); p <-  length(coef(fit))
  data_res <- data.frame(r_std = rstandard(fit), r_std_d = rstudent(fit))
  out_std <- which(abs(data_res$r_std) > t_critical)
  out_std_d <- which(abs(data_res$r_std_d) > qt(1 - alpha/(2), df = n - p - 1))
  out_h <- which(hatvalues(fit) > 2*p/n) 
  return(list('out_std' = out_std, 'out_std_d' = out_std_d, 'out_h' = out_h))
}

ponto_influente <- function(fit){
  n <- length(residuals(fit)); p <- length(coef(fit))
  i_dffits <- which(dffits(fit) >= 2*sqrt(p/n))
  i_dfbetas <- vector("list", length = p)
  for(j in 1:p){i_dfbetas[[j]] <- which(dfbetas(fit)[,j] >= 2/sqrt(n))}
  i_dist_cook = which(cooks.distance(fit) > qf(p = 0.5, df1 = p, df2 = n - p))
  return(list('i_dffits' = i_dffits, 'i_dfbetas' = i_dfbetas, 'i_dist_cook' = i_dist_cook))
}

plot(ajuste$residuals)
out <- outlier(ajuste, t_critical = 2)
pi <- ponto_influente(ajuste)

out_paises <- data.frame(a = c(paste(paises[out$out_std][1], 
                                     paises[out$out_std][2], 
                                     paises[out$out_std][3], sep = ', '), 
                         paste(paises[out$out_std_d][1],
                               paises[out$out_std_d][2],
                               paises[out$out_std_d][3], sep = ', '), 
                         paste(paises[out$out_h][1], 
                               paises[out$out_h][2], sep = ', ')))
row.names(out_paises) = c('\textit{Student}', 
                          '\textit{Student} Deletado', 'Valores Chapéu')
xtable::xtable(out_paises)

influentes <- data.frame(a = c(paste(paises[pi$i_dffits][1],
                                     paises[pi$i_dffits][2], sep = ', '), 
                               paste(paises[pi$i_dfbetas[[1]]][1], sep = ', '), 
                               paste(paises[pi$i_dfbetas[[2]]][1], sep = ', '),
                               paste(paises[pi$i_dfbetas[[4]]][1], sep = ', ')))
row.names(influentes) <- c('DFFITS', 'DFBETAS $\beta_{0}$',
                           'DFBETAS $\beta_{1}$', 'DFBETAS $\beta_{4}$')
xtable::xtable(influentes)

# Estimadores robustos

L1 <-lmrob.lar(y = lcs$sr, x = lcs[,-c(1)])
M <- rlm(y = lcs$sr, x = lcs[,-c(1)], psi = psi.huber)
S <- lmrob.S(x = lcs[,-c(1)], y = lcs$sr, control = lmrob.control(psi = "bisquare"))

summary_l1 <- summary(L1)
summary_m <- summary(M)
summary_s <- summary(S)
