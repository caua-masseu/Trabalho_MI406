library(robustbase)
library(MASS)

outlier <- function(fit, t_critical = 3, alpha = 0.05){
  n <- length(residuals(fit)); p <-  length(coef(fit))
  data_res <- data.frame(r_std = rstandard(fit), r_std_d = rstudent(fit))
  out_std <- which(abs(data_res$r_std) > t_critical)
  out_std_d <- which(abs(data_res$r_std_d) > qt(1 - alpha/(2), df = n - p - 1))
  out_h <- which(hatvalues(fit) > 2*p/n) 
  return(list('out_std' = out_std, 'out_std_d' = out_std_d, 'out_h' = out_h))
}

ponto_influente <- function(fit){
  n <- length(residuals(fit)); p <- length(coef(fit))
  i_dffits <- which(dffits(fit) > 2*sqrt(p/n))
  i_dfbetas <- vector("list", length = p)
  for(j in 1:p){i_dfbetas[[j]] <- which(dfbetas(fit)[,j] > 2/sqrt(n))}
  i_dist_cook = which(cooks.distance(fit) > qf(p = 0.5, df1 = p, df2 = n - p))
  return(list('i_dffits' = i_dffits, 'i_dfbetas' = i_dfbetas, 'i_dist_cook' = i_dist_cook))
}

#Simulacao
n_outliers_y <- 3; n_outliers_x <- 1; n_outliers <- n_outliers_x + n_outliers_y
size = 10 
size_x = 5
n <- 100 
obs_out_y <- 1:n_outliers_y
obs_out_x <- n - 1:n_outliers_x
obs_out <- c(obs_out_y, obs_out_x)

'%!in%' <- function(x,y)!('%in%'(x,y))

n_rep <- 1000
bhatmqo <- bhatl1 <- bhatm <- bhats <- data.frame()
per_res_std_y <- per_res_std_d_y <- per_res_std_x <- per_res_std_d_x <- per_h_ii <- per_res_std_false <- per_res_std_d_false <- per_h_ii_false <- c() 
per_i_dffits <- per_i_dfbeta_1 <- per_i_dfbeta_2 <- per_i_dfbeta_3 <- per_i_cook <- c()
for(i in 1:n_rep){
  set.seed(i)
  # parâmetros
  beta <- c(1, 2.5, -1.5)
  # Simulando Modelo
  x1 <- rnorm(n = n, sd = 1, mean = 3)
  x2 <- rnorm(n = n, sd = 1, mean = 4)
  X <- matrix(c(rep(1, n), x1, x2), nrow = n)
  y <- X %*% beta + rnorm(n, sd = 1)
  
  # outliers em Y
  y[obs_out] <- size + y[obs_out]
  
  #outlier em X
  x1[obs_out_x] <- x1[obs_out_x] + size_x
  
  # estimadores
  MQO <- lm(y ~ 1 + x1 + x2)
  
  #Detecao de outliers
  det_obs_out <- outlier(fit = MQO)
  
  per_res_std_false <- c(per_res_std_false, sum(det_obs_out$out_std %!in% obs_out) / (n - n_outliers) * 100)
  per_res_std_d_false <- c(per_res_std_d_false, sum(det_obs_out$out_std %!in% obs_out) / (n - n_outliers) * 100)
  per_res_std_y <- c(per_res_std_y, sum(det_obs_out$out_std %in% obs_out_y) / n_outliers_y * 100)
  per_res_std_d_y <- c(per_res_std_d_y, sum(det_obs_out$out_std_d %in% obs_out_y) / n_outliers_y * 100)
  per_res_std_x <- c(per_res_std_x, sum(det_obs_out$out_std %in% obs_out_x) / n_outliers_x * 100)
  per_res_std_d_x <- c(per_res_std_d_x, sum(det_obs_out$out_std_d %in% obs_out_x) / n_outliers_x * 100)
  per_h_ii <- c(per_h_ii, sum(det_obs_out$out_h %in% obs_out_x)/n_outliers_x * 100)
  per_h_ii_false <- c(per_h_ii_false, sum(det_obs_out$out_h %!in% obs_out_x)/(n- n_outliers_x) * 100)
  
  
  #Ponto Influente
  det_point_inf <- ponto_influente(MQO)
  per_i_dffits <- c(per_i_dffits, sum(det_point_inf$i_dffits %in% obs_out) / n_outliers * 100)
  per_i_dfbeta_1 <- c(per_i_dfbeta_1, sum(det_point_inf$i_dfbetas[[1]] %in% obs_out) / n_outliers * 100)
  per_i_dfbeta_2 <- c(per_i_dfbeta_2, sum(det_point_inf$i_dfbetas[[2]] %in% obs_out) / n_outliers * 100)
  per_i_dfbeta_3 <- c(per_i_dfbeta_3, sum(det_point_inf$i_dfbetas[[3]] %in% obs_out) / n_outliers * 100)
  per_i_cook <- c(per_i_cook, sum(det_point_inf$i_dist_cook %in% obs_out) / n_outliers * 100)
  
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

data_out <- data.frame(std = c(mean(per_res_std_false), mean(per_res_std_y), mean(per_res_std_x)), 
           std_d = c(mean(per_res_std_d_false), mean(per_res_std_d_y), mean(per_res_std_d_x)),
          h_ii = c(mean(per_h_ii_false), mean(per_h_ii), "NA"))

xtable::xtable(data_out)

data_point_influente <- data.frame(dffits = mean(per_i_dffits), dfbeta_1 = mean(per_i_dfbeta_1), dfbeta_2 = mean(per_i_dfbeta_2),
                                   dfbeta_3 = mean(per_i_dfbeta_3), cook = mean(per_i_cook))

xtable::xtable(data_point_influente)
xtable::xtable(vies)

plot(y ~ x1)
hist(y)

lm(y ~ 1 + x1 + x2)
lmrob.lar(y = y, x = X)
MASS::rlm(y = y, x = X)
lmrob.S(x = X, y = y, control = lmrob.control("KS2011", max.it = 1000))[1]
lmrob(y ~ X - 1)



