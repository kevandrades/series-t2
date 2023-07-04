################
# TP2 - SERIES #
################

##############
# INTRODUCAO #
##############

# Pacotes #
library(pacman)
p_load(tidyverse,tseries,forecast,Mcomp,kableExtra)

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
  plot(main="S?rie temporal do ?ndice geral da produ??o industrial,\n1958-1965 ",xlab="Ano",ylab="?ndice")

# Decomposicao #
  # STl
stl(serie,s.window = 12) %>% 
  plot(main= 'janela = 12') #s.window: tamanho da janela de ajuste para sazonalidade

stl(serie,s.window = 'periodic') %>% 
  plot(main= 'janela = periodic')

  # MSTL
mstl(serie, lambda = 'auto') %>% 
  autoplot() + 
  labs(x="Ano")+
  theme_bw()

#################
# MODELOS ARIMA #
#################

# Estacionariedade #
kpss.test(serie)

### MODELO            ####
### SEM TRANSFORMACAO ####

serie %>% ndiffs() 
diff(serie) %>% nsdiffs()
serie_diff <- serie %>% 
  diff() %>% 
  diff(lag=12)

kpss.test(serie_diff) 
# teste aponta estacionariedade

# ACF e PACF #
ggtsdisplay(serie_diff)
# lags simples
# ACF 
acf(serie_diff,lag.max = 12*5)
  # as correlacoes nao parecem decair para zero de forma amortizada
  # mas tbm nao parece haver quebra no lag simples  
# PACF
pacf(serie_diff, lag.max = 12*5)
  # o padrao com relacao ao lag simples nao eh mto claro 
  # o decaimento nao eh amortizado
  # nao ha quebra abrupta para o zero
  # uma vez que nao eh amortizado, testar para p = 0 e 1, P = 0 e 1; 

# conclusao: testar valores de p e q para o ARIMA

# lags sazonais
# ACF 
  # parece haver correlacao significativa no lag sazonal 1...
  # seguida de quebra 
# PACF
  # nao eh poss?vel observar decaimento amortizado
  # P=0 e Q=1 (nao identificada quebra no lag P) eh um bom candidato
  # Tbm nao identificamos quebras nos lags sazonais, pois todos estao abaixo dos limites
  # Mas nao faz sentido usar P=0 e Q=0 pq no ACF tem quebra no lag sazonal 1
  # P=1 e Q=1? nao verificado decaimento amortizado em ACF e PACF

# Testando combinacoes de p,q

melhor_AICc = Inf
for(p in 0:3) {
  for(q in 0:3) {
    fit = Arima(serie_diff, order=c(p,1,q), seasonal=c(0,1,1))
    if(fit$aicc < melhor_AICc) {
      melhor_AICc <- fit$aicc
    }
    cat("p =",p,", q =",q,", AICc =", fit$aicc, "\n")
  }
}
melhor_AICc # p=0 q=2 AICc = 881.2515


# Ajuste do modelo escolhido
(fit = Arima(serie_diff, order=c(0,1,2), seasonal=c(0,1,1)))

# Residuos
par(mfrow=c(1,1))
residuos <- fit$residuals %>% window(start=c(1960,1))
par(mfrow=c(1,3))
plot(residuos,main="Res?duos ap?s \ninicializa??o do modelo")
qqnorm(residuos);qqline(residuos)
acf(residuos,lag.max=12*5)
shapiro <- shapiro.test(residuos)
kpss <- kpss.test(residuos)
box <- Box.test(residuos,lag=15,type='Ljung-Box')

# resumo dos testes
df.residuos <- data.frame("Teste"=c("Shapiro-Wilk","KPSS","Ljung-Box"),
                        "P-valor"=c(shapiro$p.value,kpss$p.value,box$p.value),
                        check.names=FALSE)



### MODELO            ####
### COM TRANSFORMACAO ####
### DE BOX-COX        ###

lambda <- BoxCox.lambda(serie)
paste("Lambda =", lambda)

serie_transf <- BoxCox(serie,lambda = lambda)
par(mfrow=c(1,2))
plot.ts(serie,main="S?rie Original")
plot.ts(serie_transf,main="S?rie Transformada")

serie_transf %>% ndiffs()
serie_transf %>% diff() %>% nsdiffs()

serie_transf_diff <- serie_transf %>% 
  diff() %>% diff(lag=12)

kpss.test(serie_transf_diff)

# ACF e PACF #
par(mfrow=c(1,2))
acf(serie_transf_diff,lag.max = 5*12)
pacf(serie_transf_diff,lag.max = 5*12)


#################
# MODELOS ETS #
#################

library(tseries)
library(tidyverse)
library(forecast)
library(knitr)

# ETS
# Resultado de critério de informação ETS sem transformação
fit1<- ets(serie,model = "AAA")
fit2<- ets(serie,model = "AAA",damped = TRUE)
fit3<- ets(serie,model = "MAA")
fit4<- ets(serie,model = "MAA",damped = TRUE)
fit5<- ets(serie,model = "MAM")
fit6<- ets(serie,model = "MMM")
fit7<- ets(serie,model = "MAM",damped = TRUE)
fit8<- ets(serie,model = "MMM", damped = TRUE)
AIC <- rbind(fit1$aic,fit2$aic,fit3$aic,fit4$aic,
             fit5$aic,fit6$aic,fit7$aic,fit8$aic)
AICc <- rbind(fit1$aicc,fit2$aicc,fit3$aicc,fit4$aicc,
              fit5$aicc,fit6$aicc,fit7$aicc,fit8$aicc)
BIC <- rbind(fit1$bic,fit2$bic,fit3$bic,fit4$bic,
             fit5$bic,fit6$bic,fit7$bic,fit8$bic)
Modelo <- cbind(c("ETS(A,A,A)","ETS(A,Ad,A)","ETS(M,A,A)","ETS(M,Ad,A)",
                  "ETS(M,A,M)","ETS(M,M,M)","ETS(M,Ad,M)","ETS(M,Md,M)"))
d <- data.frame(Modelo,AIC,AICc,BIC)
knitr::kable(d) # MODELO ETS(M,M,M) MENOR AIC, AICc e BIC


# Decomposição ETS sem transformação
plot(fit6)

# Análise de resíduos ETS sem transformação
E <- fit6$residuals
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
