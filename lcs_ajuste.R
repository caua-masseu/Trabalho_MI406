data("LifeCycleSavings")

library(ggplot2)
library(tidyverse)
library(GGally)

library(robustbase)
library(MASS)

# Analise exploratoria

ggpairs(LifeCycleSavings[,])
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

# Teste Outlier
out <- outlier(ajuste, t_critical = 2)

# Teste Ponto Influente
pi <- ponto_influente(ajuste)

# Tabelas Detectados
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

# Detectados
outliers <- unique(c(out$out_std, out$out_std_d, out$out_h))
pi <- unique(c(pi$i_dffits, pi$i_dist_cook, pi$i_dfbetas[[1]], pi$i_dfbetas[[2]],
         pi$i_dfbetas[[4]]))
intersect(outliers, pi)

# Exemplos de diagnosticos
# plot(ajuste)

# Residuos e ajustados do MQO
ajuste_data <- data.frame(fitted = ajuste$fitted.values, 
                          residuos = ajuste$residuals, ajuste$model,
                          pais = row.names(lcs))
ajuste_data$outlier <- 'Nada Detectado'
ajuste_data$outlier[outliers] <- 'Outlier'
ajuste_data$outlier[pi] <- 'Outlier e Ponto Influente'
ajuste_data$outlier <- factor(ajuste_data$outlier, 
                                 levels  = c('Outlier e Ponto Influente', 
                                             'Outlier', 'Nada Detectado'))

# Residuos x Ajustados
ggplot() + 
  geom_point(data = ajuste_data, aes(x = fitted, y = residuos, color = outlier), size = 5) + 
  # geom_point(data = ajuste_data[c('Japan', 'Libya', 'Zambia'),], 
  #            aes(x = fitted, y = residuos), color = 'red') + 
  # geom_point(data = ajuste_data[c('Chile', 'United States'),], 
  #            aes(x = fitted, y = residuos), color = 'blue') + 
  geom_text(data = ajuste_data[c('Australia', 'Chile', 'Japan', 'Libya',  
                                 'Malaysia', 'United States',
                                 'Zambia'),], 
            aes(x = fitted, y = residuos), size = 5, 
            label = c('Austrália', 'Chile', 'Japão', 'Líbia', 'Malásia', 'EUA', 'Zâmbia'),nudge_x = -0.4) + 
  theme_bw(base_size = 20) + 
  scale_color_manual(values = c('red', 'blue', 'black')) + 
  labs(y = 'Resíduos', x = 'Valores Ajustados', color = 'Detecção', 
       subtitle = 'Estimador: MQO') + 
  ylim(c(-8, 11)) + xlim(c(5, 14))

# Reta Residuo
ggplot() +
  geom_point(data = ajuste_data, aes(x = fitted, y = sr, 
                                 color = outlier), size = 5) +
  geom_text(data = ajuste_data[c('Australia', 'Chile', 'Japan', 'Libya', 'Malaysia', 'United States',
                             'Zambia'),],
            aes(x =  fitted, y = sr), size = 5,
            label = c('Austrália', 'Chile', 'Japão', 'Líbia', 'Malásia', 'EUA', 'Zâmbia'),nudge_x = -0.4) +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c('red', 'blue', 'black')) +
  labs(y = 'Resposta', x = 'Valores Ajustados', color = 'Detecção', 
       subtitle = 'Estimador: MQO') +
  geom_line(data = ajuste_data, aes(x = sr, y = sr), linewidth = 1.5, linetype = 'dashed')


# QQPlot
g <- ggplot(data = ajuste_data, aes(sample = residuos)) +
  geom_qq(aes(sample = residuos)) + 
  geom_qq_line()
qq <- ggplot_build(g)
qq_coord <- data.frame(qq$data[[1]][,c(1,2)], qq$plot$data)

# Comparacao entre ggplot e ggplot reconstruido para colorir
# plot(g)
# ggplot() +
#   geom_qq_line(data = qq_coord, aes(sample = residuos)) +
#   geom_point(data = qq_coord, aes(x = x, y = y, color = outlier)) +
#   theme_bw() +
#   geom_text(data = qq_coord[c('Chile', 'Japan', 'Libya', 'United States',
#                                  'Zambia'),],
#             aes(x = x, y = y),
#             label = c('Chile', 'Japão', 'Líbia', 'EUA', 'Zâmbia'),nudge_x = -0.3,
#             nudge_y = 0.4) +
#   scale_color_manual(values = c('red', 'blue', 'black')) +
#   labs(color = 'Detecção') +
#   theme(legend.position = 'none')

ggplot() +
  geom_qq_line(data = qq_coord, aes(sample = residuos)) +
  geom_point(data = qq_coord, aes(x = x, y = y, color = outlier), size = 3) +
  theme_bw(base_size = 20) +
  geom_text(data = qq_coord[c('Australia', 'Chile', 'Japan', 'Libya', 'Malaysia', 'United States',
                                 'Zambia'),],
            aes(x = x, y = y), size = 5,
            label = c('Austrália', 'Chile', 'Japão', 'Líbia', 'Malásia', 
                      'EUA', 'Zâmbia'),nudge_x = -0.3,
            nudge_y = 0.4) +
  scale_color_manual(values = c('red', 'blue', 'black')) +
  labs(color = 'Detecção', x = 'Quantis Teóricos', y = 'Quantis Amostrais', 
       subtitle = 'Estimador: MQO') + 
  ylim(c(-8.75, 12)) + xlim(c(-2.75, 2.75))

# Estimadores robustos
X_intercepto <- cbind(intercepto = c(rep(1, 50)), lcs[,-c(1)])
L1 <-lmrob.lar(y = lcs$sr, x = X_intercepto)
M <- rlm(y = lcs$sr, x = X_intercepto, psi = psi.huber)
S <- lmrob.S(x = X_intercepto, y = lcs$sr, control = lmrob.control(psi = "bisquare"))

# L1
L1_data <- data.frame(cbind(L1$residuals, fitted = (lcs$sr - L1$residuals), 
                            ajuste_data$pais,lcs, ajuste_data$outlier))
# Residuos x Ajustados
ggplot() + 
  geom_point(data = L1_data, aes(x = fitted, y = L1.residuals, 
                                 color = ajuste_data.outlier), 
             size = 5) + 
  geom_text(data = L1_data[c('Australia', 'Chile', 'Japan', 'Libya', 'Malaysia', 'United States',
                                 'Zambia'),], 
            aes(x = fitted, y = L1.residuals), size = 5, 
            label = c('Austrália', 'Chile', 'Japão', 'Líbia','Malásia', 'EUA', 'Zâmbia'),nudge_x = -0.4) + 
  theme_bw(base_size = 20) + 
  scale_color_manual(values = c('red', 'blue', 'black')) + 
  labs(y = 'Resíduos', x = 'Valores Ajustados', color = 'Detecção', 
       subtitle = 'Estimador: L1') + 
  ylim(c(-8, 14)) + xlim(c(3, 14))

# Reta Residuo
ggplot() +
  geom_point(data = L1_data, aes(x = fitted, y = sr, 
                                 color = ajuste_data.outlier), size = 5) +
  geom_text(data = L1_data[c('Australia', 'Chile', 'Japan', 'Libya', 'Malaysia', 'United States',
                            'Zambia'),],
            aes(x =  fitted, y = sr), size = 5,
            label = c('Austrália', 'Chile', 'Japão', 'Líbia', 'Malásia', 'EUA', 'Zâmbia'),nudge_x = -0.4) +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c('red', 'blue', 'black')) +
  labs(y = 'Resposta', x = 'Valores Ajustados', color = 'Detecção', 
       subtitle = 'Estimador: L1') +
  geom_line(data = L1_data, aes(x = sr, y = sr), linewidth = 1.5, linetype = 'dashed')

# QQPlot
g <- ggplot(data = L1_data, aes(sample = L1.residuals)) +
  geom_qq(aes(sample = L1.residuals)) + 
  geom_qq_line()
qq <- ggplot_build(g)
qq_coord <- data.frame(qq$data[[1]][,c(1,2)], qq$plot$data)

ggplot() +
  geom_qq_line(data = qq_coord, aes(sample = L1.residuals)) +
  geom_point(data = qq_coord, aes(x = x, y = y, color = ajuste_data.outlier), size = 3) +
  theme_bw(base_size = 20) +
  geom_text(data = qq_coord[c('Australia', 'Chile', 'Japan', 'Libya', 'Malaysia', 'United States',
                              'Zambia'),],
            aes(x = x, y = y), size = 5,
            label = c('Austrália', 'Chile', 'Japão', 'Líbia', 'Malásia', 
                      'EUA', 'Zâmbia'),nudge_x = -0.3,
            nudge_y = 0.4) +
  scale_color_manual(values = c('red', 'blue', 'black')) +
  labs(color = 'Detecção', x = 'Quantis Teóricos', y = 'Quantis Amostrais', 
       subtitle = 'Estimador: L1') + 
  ylim(c(-8.75, 13.5)) + xlim(c(-2.75, 2.75))

# M
M_data <- data.frame(cbind(M$residuals, M$fitted.values, lcs, ajuste_data$pais, outlier = ajuste_data$outlier))
# Residuos x Ajustados
ggplot() + 
  geom_point(data = M_data, aes(x = M.fitted.values, y = M.residuals, color = outlier), size = 5) + 
  geom_text(data = M_data[c('Australia','Chile', 'Japan', 'Libya','Malaysia', 'United States',
                                 'Zambia'),], 
            aes(x =  M.fitted.values, y = M.residuals), size = 5, 
            label = c('Austrália', 'Chile', 'Japão', 'Líbia','Malásia', 'EUA', 'Zâmbia'),nudge_x = -0.4) + 
  theme_bw(base_size = 20) + 
  scale_color_manual(values = c('red', 'blue', 'black')) + 
  labs(y = 'Resíduos', x = 'Valores Ajustados', color = 'Detecção', 
       subtitle = 'Estimador: M')+
  ylim(c(-8, 14)) + xlim(c(3, 14))

# Reta residuo
ggplot() +
  geom_point(data = M_data, aes(x = M.fitted.values, y = sr, color = outlier), size = 5) +
  geom_text(data = M_data[c('Australia', 'Chile', 'Japan', 'Libya', 'Malaysia', 'United States',
                            'Zambia'),],
            aes(x =  M.fitted.values, y = sr), size = 5,
            label = c('Austrália', 'Chile', 'Japão', 'Líbia', 'Malásia', 'EUA', 'Zâmbia'),nudge_x = -0.4) +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c('red', 'blue', 'black')) +
  labs(y = 'Resposta', x = 'Valores Ajustados', color = 'Detecção', 
       subtitle = 'Estimador: M') +
  geom_line(data = M_data, aes(x = sr, y = sr), linewidth = 1.5, linetype = 'dashed')

# QQPlot
g <- ggplot(data = M_data, aes(sample = M.residuals)) +
  geom_qq(aes(sample = M.residuals)) + 
  geom_qq_line()
qq <- ggplot_build(g)
qq_coord <- data.frame(qq$data[[1]][,c(1,2)], qq$plot$data)

# Reta residuo
ggplot() +
  geom_point(data = M_data, aes(x = M.fitted.values, y = sr, color = M$w), size = 5) +
  geom_text(data = M_data[c('Australia', 'Chile', 'Japan', 'Libya', 'Malaysia', 'United States',
                            'Zambia'),],
            aes(x =  M.fitted.values, y = sr), size = 5,
            label = c('Austrália', 'Chile', 'Japão', 'Líbia', 'Malásia', 'EUA', 'Zâmbia'),nudge_x = -0.4) +
  theme_bw(base_size = 20) +
  labs(y = 'Resposta', x = 'Valores Ajustados', color = 'Pesos em MQPI', 
       subtitle = 'Estimador: M') +
  geom_line(data = M_data, aes(x = sr, y = sr), linewidth = 1.5, linetype = 'dashed')

ggplot() +
  geom_qq_line(data = qq_coord, aes(sample = M.residuals)) +
  geom_point(data = qq_coord, aes(x = x, y = y, color = outlier), size = 3) +
  theme_bw(base_size = 20) +
  geom_text(data = qq_coord[c('Australia', 'Chile', 'Japan', 'Libya', 'Malaysia', 'United States',
                              'Zambia'),],
            aes(x = x, y = y), size = 5,
            label = c('Austrália', 'Chile', 'Japão', 'Líbia', 'Malásia', 
                      'EUA', 'Zâmbia'),nudge_x = -0.3,
            nudge_y = 0.4) +
  scale_color_manual(values = c('red', 'blue', 'black')) +
  labs(color = 'Detecção', x = 'Quantis Teóricos', y = 'Quantis Amostrais', 
       subtitle = 'Estimador: M') + 
  ylim(c(-8.75, 13.5)) + xlim(c(-2.75, 2.75))


# S
S_data <- data.frame(cbind(S$residuals, S$fitted.values, lcs, ajuste_data$pais, outlier = ajuste_data$outlier))
# Residuos x Ajustados
ggplot() +
  geom_point(data = S_data, aes(x = S.fitted.values, y = S.residuals, color = outlier), size = 5) +
  # geom_point(data = ajuste_data[c('Japan', 'Libya', 'Zambia'),],
  #            aes(x = fitted, y = residuos), color = 'red') +
  # geom_point(data = ajuste_data[c('Chile', 'United States'),],
  #            aes(x = fitted, y = residuos), color = 'blue') +
  geom_text(data = S_data[c('Australia','Chile', 'Japan', 'Libya', 'Malaysia', 'United States',
                            'Zambia'),],
            aes(x =  S.fitted.values, y = S.residuals), size = 5,
            label = c('Austrália', 'Chile', 'Japão', 'Líbia', 'Malásia', 'EUA', 'Zâmbia'),nudge_x = -0.4) +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c('red', 'blue', 'black')) +
  labs(y = 'Resíduos', x = 'Valores Ajustados', color = 'Detecção', 
       subtitle = 'Estimador: S') +
  ylim(c(-8, 14)) + xlim(c(3, 14))

# Reta residuo
ggplot() +
  geom_point(data = S_data, aes(x = S.fitted.values, y = sr, color = outlier), size = 5) +
  geom_text(data = S_data[c('Australia', 'Chile', 'Japan', 'Libya', 'Malaysia', 'United States',
                            'Zambia'),],
            aes(x =  S.fitted.values, y = sr), size = 5,
            label = c('Austrália','Chile', 'Japão', 'Líbia', 'Malásia', 'EUA', 'Zâmbia'),nudge_x = -0.4) +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c('red', 'blue', 'black')) +
  labs(y = 'Resposta', x = 'Valores Ajustados', color = 'Detecção', 
       subtitle = 'Estimador: S') +
  geom_line(data = S_data, aes(x = sr, y = sr), linewidth = 1.5, linetype = 'dashed')


# QQPlot
g <- ggplot(data = S_data, aes(sample = S.residuals)) +
  geom_qq(aes(sample = S.residuals)) + 
  geom_qq_line()
qq <- ggplot_build(g)
qq_coord <- data.frame(qq$data[[1]][,c(1,2)], qq$plot$data)

ggplot() +
  geom_qq_line(data = qq_coord, aes(sample = S.residuals)) +
  geom_point(data = qq_coord, aes(x = x, y = y, color = outlier), size = 3) +
  theme_bw(base_size = 20) +
  geom_text(data = qq_coord[c('Australia', 'Chile', 'Japan', 'Libya', 'Malaysia', 'United States',
                              'Zambia'),],
            aes(x = x, y = y), size = 5,
            label = c('Austrália', 'Chile', 'Japão', 'Líbia', 'Malásia', 
                      'EUA', 'Zâmbia'),nudge_x = -0.3,
            nudge_y = 0.4) +
  scale_color_manual(values = c('red', 'blue', 'black')) +
  labs(color = 'Detecção', x = 'Quantis Teóricos', y = 'Quantis Amostrais', 
       subtitle = 'Estimador: S') + 
  ylim(c(-8.75, 13.5)) + xlim(c(-2.75, 2.75))

# Reta residuo
ggplot() +
  geom_point(data = S_data, aes(x = S.fitted.values, y = sr, color = S$rweights), size = 5) +
  geom_text(data = S_data[c('Australia', 'Chile', 'Japan', 'Libya', 'Malaysia', 'United States',
                            'Zambia'),],
            aes(x =  S.fitted.values, y = sr), size = 5,
            label = c('Austrália', 'Chile', 'Japão', 'Líbia', 'Malásia', 'EUA', 'Zâmbia'),nudge_x = -0.4) +
  theme_bw(base_size = 20) +
  labs(y = 'Resposta', x = 'Valores Ajustados', color = 'Pesos em MQPI', 
       subtitle = 'Estimador: S') +
  geom_line(data = S_data, aes(x = sr, y = sr), linewidth = 1.5, linetype = 'dashed')



# Coeficientes
ajuste$coefficients
L1$coefficients
M$coefficients
S$coefficients
(M$w)
which(S$rweights == 0)
S$scale
