################
# TP2 - SERIES #
################

# ----------------------------------------------------------- #
# 1. INTRODUCAO #
# ----------------------------------------------------------- #

# Pacotes #
library(pacman)
p_load(knitr,tidyverse,tseries,forecast,Mcomp,kableExtra, purrr)

# Serie escolhida #
data(M3)
id=2507

serie <- M3[[id]]$x # dados de treinamento
dados_teste <- M3[[id]]$xx # dados de teste
horizonte <- M3[[id]]$h

M3[[id]] # infos gerais serie
M3[[id]]$description # descricao
M3[[id]]$n # no. observacoes
M3[[id]]$period # periodicidade
M3[[id]]$type #tipo

serie %>%
  plot(
    main = "Série temporal do ?ndice geral da produ??o industrial,\n1958-1965 ",
    xlab = "Ano",
    ylab = "Índice"
  )

# Decomposicao #
# STl
stl(serie, s.window = 12) %>%
  plot(main = "janela = 12")
#s.window: tamanho da janela de ajuste para sazonalidade

stl(serie, s.window = "periodic") %>%
  plot(main = "janela = periodic")

# MSTL
mstl(serie, lambda = "auto") %>%
  autoplot() +
  labs(x = "Ano") +
  theme_bw()

# ----------------------------------------------------------- #
# 2. MODELOS ARIMA
# ----------------------------------------------------------- #
# Estacionariedade #
kpss.test(serie)

### 2.1 Modelo sem transformação

serie %>% ndiffs()
diff(serie) %>% nsdiffs()
serie_diff <- serie %>%
  diff() %>%
  diff(lag = 12)

kpss.test(serie_diff) 
# teste aponta estacionariedade

# ACF e PACF #
ggtsdisplay(serie_diff)
# lags simples
# ACF 
acf(serie_diff, lag.max = 12 * 5)
# PACF
pacf(serie_diff, lag.max = 12*5)

# conclusao: testar valores de p e q para o ARIMA

# Testando combinacoes de p,q

melhor_AICc <- Inf

for(p in 0:3) {
  for(q in 0:3) {
    fit <- Arima(serie_diff, order=c(p, 0, q), seasonal=c(0, 0, 1))
    if(fit$aicc < melhor_AICc) {
      melhor_AICc <- fit$aicc
    }
    cat("p =", p,", q =", q,", AICc =", fit$aicc, "\n")
  }
}
melhor_AICc # p=0 q=1 AICc = 1001.2

# Ajuste do modelo escolhido
(fit_arima <- Arima(serie_diff, order=c(0, 0, 1), seasonal = c(0, 0, 1)))

fit_arima$aicc 

# Residuos
par(mfrow = c(1, 1))
residuos <- fit_arima$residuals %>%
  window(start = c(1960, 1))
par(mfrow=c(1, 3))
plot(residuos, main = "Resíduos após ininicialização do modelo")
qqnorm(residuos)
qqline(residuos)
acf(residuos, lag.max = 12 * 5)
shapiro <- shapiro.test(residuos)
kpss <- kpss.test(residuos)
box <- Box.test(residuos,lag=15,type="Ljung-Box")

# resumo dos testes
df.residuos <- data.frame(
  Teste = c("Shapiro-Wilk", "KPSS", "Ljung-Box"),
  "P-valor" = c(shapiro$p.value, kpss$p.value, box$p.value),
  check.names = FALSE
)

### 1.2 Modelo com transformação de BoxCox

lambda <- BoxCox.lambda(serie)
paste("Lambda =", lambda)

serie_transf <- BoxCox(serie, lambda = lambda)
par(mfrow = c(1, 2))
plot.ts(serie, main = "Série Original")
plot.ts(serie_transf, main = "Série Transformada")

serie_transf %>% ndiffs()
serie_transf %>% diff() %>% nsdiffs()

serie_transf_diff <- serie_transf %>%
  diff() %>%
  diff(lag=12)

kpss.test(serie_transf_diff)

# ACF e PACF #
par(mfrow=c(1,2))
acf(serie_transf_diff,lag.max = 5*12)
pacf(serie_transf_diff,lag.max = 5*12)

## Testar combinacoes de p,q,P,Q ##

melhor_AICc <- Inf
melhores_qs <- melhores_ps <- melhores_Ps <- melhores_Qs <- AICcs <- c()

for (q in 0:3) {
  for (p in 0:3) {
    for (P in 0:1) {
      for (Q in 0:1) {
        deu_errado <- FALSE
        
        tryCatch({
          fit2 <- Arima(serie_transf_diff, order = c(p, 0, q), seasonal = c(P, 0, Q))
          melhor_AICc <- fit2$aicc
          melhores_ps <- c(melhores_ps, p)
          melhores_qs <- c(melhores_qs, q)
          melhores_Ps <- c(melhores_Ps, P)
          melhores_Qs <- c(melhores_Qs, Q)
          AICcs <- c(AICcs, melhor_AICc)
        }, error = function(e) {
          deu_errado <- TRUE
        })
        
        if (deu_errado) {
          next
        }
      }
    }
  }
}

aiccs<- data.frame(
  "AICc"=AICcs,
  "p"=melhores_ps,
  "q"=melhores_qs,
  "P"=melhores_Ps,
  "Q"=melhores_Qs
)

aiccs %>%
  arrange(desc(AICc)) %>%
  slice(59:64) %>%
  kableExtra::kbl(.,digits=4,align=c("l","c","c"),booktabs = T,caption = "AICc dos modelos") %>% 
  kableExtra::kable_classic(full_width=FALSE,latex_options = "HOLD_position")

# menor AICc = -539.9822 (0,0,1) x (0,0,1)

# Ajuste do modelo escolhido
(fit_arima_boxcox <- Arima(
  serie_transf_diff, order=c(0, 0, 1), seasonal=c(0, 0, 1))
)

fit_arima_boxcox$aicc 

# Residuos
par(mfrow=c(1, 1))
residuos <- fit_arima_boxcox$residuals %>% window(start=c(1960,1))
par(mfrow=c(1, 3))
plot(residuos, main="Resíduos após ininicialização do modelo")
qqnorm(residuos)
qqline(residuos)
acf(residuos,lag.max = 12 * 5)
shapiro <- shapiro.test(residuos)
kpss <- kpss.test(residuos)
box <- Box.test(residuos, lag = 15, type="Ljung-Box")

# resumo dos testes
df.residuos <- data.frame(
  Teste = c("Shapiro-Wilk", "KPSS", "Ljung-Box"),
  "P-valor" = c(shapiro$p.value, kpss$p.value, box$p.value),
  check.names=FALSE
)
#################
# MODELOS ETS #
#################

# ETS
# Resultado de critério de informação ETS sem transformação
ets_fit1 <- ets(serie,model = "AAA")
ets_fit2 <- ets(serie,model = "AAA", damped = TRUE)
ets_fit3 <- ets(serie,model = "MAA")
ets_fit4 <- ets(serie,model = "MAA",damped = TRUE)
ets_fit5 <- ets(serie,model = "MAM")
ets_fit6 <- ets(serie,model = "MMM")
ets_fit7 <- ets(serie,model = "MAM",damped = TRUE)
ets_fit8 <- ets(serie,model = "MMM", damped = TRUE)

d <- list(ets_fit1, ets_fit2, ets_fit3, ets_fit4, ets_fit5, ets_fit6, ets_fit7, ets_fit8) %>%
  map(
    ~with(., c(Method = method, AIC = aic, AICc = aicc, BIC = bic))
  ) %>%
  reduce(rbind) %>%
  as.data.frame(row.names=FALSE)

knitr::kable(d) # MODELO ETS(M,M,M) MENOR AIC, AICc e BIC


# Decomposição ETS sem transformação
plot(ets_fit6)

# Análise de resíduos ETS sem transformação
E <- ets_fit6$residuals
par(mfrow=c(2,2))
plot(E)
acf(E)
pacf(E)
qqnorm(E)
qqline(E)

# Testes para ETS sem transformação
p_valor <- c(shapiro.test(E)$p.value,kpss.test(E)$p.value,
             Box.test(E,lag=15,type="Ljung-Box",fitdf=3)$p.value)
Estatistica <- c(shapiro.test(E)$statistic,kpss.test(E)$statistic,
                 Box.test(E,lag=15,type="Ljung-Box",fitdf=3)$statistic)
Teste <- c("Normalidade","Estacionariedade","Independencia")
d <- data.frame(Estatistica,p_valor)
knitr::kable(d)

# ETS com transformação

lambda <- serie %>% BoxCox.lambda()
serie_box <- serie %>% BoxCox(lambda)

# Visualização e decomposição da ETS com transformação
dev.off
plot(serie_box,main="Série com\ntransformacao de Box-Cox")
mstl(serie_box)%>%plot()

# Resultado de critério de informação ETS com transformação
ets_fit1b<- ets(serie_box,model = "AAA")
ets_fit2b<- ets(serie_box,model = "AAA",damped = TRUE)
ets_fit3b<- ets(serie_box,model = "MAA")
ets_fit4b<- ets(serie_box,model = "MAA",damped = TRUE)
ets_fit5b<- ets(serie_box,model = "MAM")
ets_fit6b<- ets(serie_box,model = "MMM")
ets_fit7b<- ets(serie_box,model = "MAM",damped = TRUE)
ets_fit8b<- ets(serie_box,model = "MMM", damped = TRUE)


d_boxcox <- list(ets_fit1b, ets_fit2b, ets_fit3b, ets_fit4b, ets_fit5b, ets_fit6b, ets_fit7b, ets_fit8b) %>%
  map(
    ~with(., c(Method = method, AIC = aic, AICc = aicc, BIC = bic))
  ) %>%
  reduce(rbind) %>%
  as.data.frame(row.names=FALSE)
knitr::kable(d_boxcox)

# modelo escolhido FIT1B - menos AIC, AICc e BIC
ets_fit1b

# ETS(A,A,A) - parâmetros
#alpha = 0.55 
#beta  = 1e-04 
#gamma = 1e-04 

# Decomposição ETS com transformação
plot(ets_fit1b)

# Análise de resíduos ETS com transformação
Eb <- ets_fit1$residuals
par(mfrow=c(2,2))
plot(Eb)
acf(Eb)
pacf(Eb)
qqnorm(Eb)
qqline(Eb)

# Testes para ETS com transformação
p_valor <- c(shapiro.test(Eb)$p.value,kpss.test(Eb)$p.value,
             Box.test(Eb,lag=15,type="Ljung-Box",fitdf=3)$p.value)
Estatística <- c(shapiro.test(Eb)$statistic,kpss.test(Eb)$statistic,
                 Box.test(Eb,lag=15,type="Ljung-Box",fitdf=3)$statistic)
Teste <- c("Normalidade","Estacionariedade","Independência")
db <- data.frame(Estatística,p_valor)
knitr::kable(db)

#################
#    PREVISÃO   #
#################

#Funções de previsão
f_arima <- function(y, h, ...){
  fit = Arima(y, order=c(0, 1, 1), seasonal = c(0, 1, 1))
  forecast(fit, h, ...)
}

# #Sarima com transformação 
f_arima_boxcox <- function(y, h, ...){
  fit = Arima(y, order=c(0, 1, 1), seasonal=c(0, 1, 1), lambda = lambda)
  forecast(fit, h, ...)
}

#ETS
f_ets <- function(y, h, ...){
  fit = ets(y, model="MMM")
  forecast(fit, h, ...)
}

#ETS com transformação
f_ets_boxcox <- function(y, h, ...){
  fit = ets(y, model="AAA", lambda = lambda)
  forecast(fit, h, ...)
}

# gráficos de previsão
par(mfrow=c(2, 2))
plot(
  f_arima(y=serie_diff, h=5, level = 95)
)

plot(
  f_arima_boxcox(y=serie_diff, h=5, level = 95)
)

plot(
  f_ets(y=serie, h=5, level = 95)
)

plot(
  f_ets_boxcox(y=serie, h=5, level = 95)
)

# Tamanho da série

n = length(serie)

# Erros de previsão

# Sarima

CV_arima = tsCV(
  y = serie, forecastfunction = f_arima,
  h = 5, initial = n - 14
)

#Sarima com transformação
CV_arima_boxcox = tsCV(
  y = serie, forecastfunction = f_arima_boxcox,
  h = 5, initial = n - 14
)

#ETS
CV_ets = tsCV(
  y = serie, forecastfunction = f_ets,
  h = 5, initial = n - 14
)

#ETS com transformação
CV_ets_boxcox = tsCV(
  y = serie, forecastfunction = f_ets_boxcox,
  h = 5, initial = n - 14
)


# MAE preditivo
MAEs_tab <- data.frame(
  MAE_arima = CV_arima %>% abs() %>% colMeans(na.rm=T),
  MAE_arima_boxcox = CV_arima_boxcox %>% abs() %>% colMeans(na.rm=T),
  MAE_ets = CV_ets %>% abs() %>% colMeans(na.rm=T),
  MAE_ets_boxcox = CV_ets_boxcox %>% abs() %>% colMeans(na.rm=T),
  h=1:5
) %>%
  pivot_longer(cols=c("MAE_arima", "MAE_arima_boxcox", "MAE_ets", "MAE_ets_boxcox")) %>%
  mutate(
    name = c(
      "MAE_arima" = "SARIMA",
      "MAE_arima_boxcox" = "SARIMA Box-Cox",
      "MAE_ets" = "ETS",
      "MAE_ets_boxcox" = "ETS Box-Cox"
    )[name]
  )


ggplot(MAEs_tab) +
  aes(x = h, y = value, color = name) +
  geom_line() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = c(0.87, 0.26),
    legend.title = element_blank(),
    legend.spacing.y = unit(0, "mm"), 
    panel.border = element_rect(colour = "black", fill=NA),
    axis.text = element_text(colour = 1, size = 12),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  ) +
  labs(x = "Janela", y = "Erro Médio Absoluto") +
  scale_color_viridis_d()

ggsave("mae_modelos.pdf")




#################
#      MAE      #
#################

# já deixei indicado a sugestão de código para a arima com transformação

#arima
xx.forec_arima <- f_arima(serie, horizonte)
#ets
xx.forec_ets <-f_ets(serie, horizonte)
#arima_boxcox
xx.forec_arima_boxcox <-f_arima_boxcox(serie, horizonte)
#ets_boxcox
xx.forec_ets_boxcox <-f_ets_boxcox(serie, horizonte)
#auto.arima
xx.forec_auto <- auto.arima(serie, allowdrift = FALSE) %>% forecast(horizonte)
#ses
xx.forec_ses <- ses(serie, allowdrift = FALSE) %>% forecast(horizonte)
#holt
xx.forec_holt <- holt(serie, allowdrift = FALSE) %>% forecast(horizonte)
#ets
xx.forec_ets <- ets(serie) %>% forecast(horizonte)
#stlf
xx.forec_stlf <- stlf(serie) %>% forecast(horizonte)
#bats
xx.forec_bats <- bats(serie, allowdrift=FALSE) %>% forecast(horizonte)
#tbats
xx.forec_tbats <- tbats(serie, allowdrift=FALSE) %>% forecast(horizonte)

## erro absoluto médio da previsão

MAE_arima2 <- mean(abs(dados_teste - xx.forec_arima$mean))
MAE_ets2 <- mean(abs(dados_teste - xx.forec_ets$mean))
MAE_arima_boxcox2 <- mean(abs(dados_teste - xx.forec_arima_boxcox$mean))
MAE_ets_boxcox2 <- mean(abs(dados_teste - xx.forec_ets_boxcox$mean))
MAE_auto <- mean(abs(dados_teste - xx.forec_auto$mean))
MAE_ses <- mean(abs(dados_teste - xx.forec_ses$mean))
MAE_holt <- mean(abs(dados_teste - xx.forec_holt$mean))
MAE_ets <- mean(abs(dados_teste - xx.forec_ets$mean))
MAE_stlf <- mean(abs(dados_teste - xx.forec_stlf$mean))
MAE_bats <- mean(abs(dados_teste - xx.forec_bats$mean))
MAE_tbats <- mean(abs(dados_teste - xx.forec_tbats$mean))
data_mae <- rbind(MAE_arima2,MAE_ets2, MAE_arima_boxcox2, MAE_ets_boxcox2,
                  MAE_auto,MAE_ets,MAE_holt,MAE_ets,MAE_stlf,MAE_bats,MAE_tbats)
data_mae <- as.data.frame(data_mae)
colnames(data_mae) <- "MAE"
rownames(data_mae) <- c("Arima", "ETS", "Arima Box Cox","ETS Box Cox",
                        "Auto arima", "Ses", "Holt", "Ets",
                        "Stlf", "Bats", "Tbats")

knitr::kable(data_mae)
