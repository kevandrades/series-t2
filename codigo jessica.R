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
  plot(main="Série temporal do Índice geral da produção industrial,\n1958-1965 ",xlab="Ano",ylab="Índice")

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
plot(residuos,main="Resíduos após ininicialização do modelo")
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
plot.ts(serie,main="Série Original")
plot.ts(serie_transf,main="Série Transformada")

serie_transf %>% ndiffs()
serie_transf %>% diff() %>% nsdiffs()

serie_transf_diff <- serie_transf %>% 
  diff() %>% diff(lag=12)

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
          fit2 <- Arima(serie_transf_diff, order = c(p, 1, q), seasonal = c(P, 0, Q))
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

aiccs<- data.frame("AICc"=AICcs,
                   "p"=melhores_ps,
                   "q"=melhores_qs,
                   "P"=melhores_Ps,
                   "Q"=melhores_Qs
)

aiccs %>%
  arrange(desc(AICc)) %>%
  slice(59:64) %>%
  kableExtra::kbl(.,digits=4,align=c('l','c','c'),booktabs = T,caption = 'AICc dos modelos') %>% 
  kableExtra::kable_classic(full_width=FALSE,latex_options = "HOLD_position")



# Ajuste do modelo escolhido
(fit2 = Arima(serie_transf_diff, order=c(0,1,3), seasonal=c(0,1,1),
              lambda = lambda))
fit2$aicc

# Obs.: em alguns cálculos (principalmente com P=0 e Q=0), apareceu o erro "non-stationary seasonal AR part from CSS". 
# Isso acontece porque, ao usar CSS (soma condicional de quadrados), é possível que os coeficientes autoregressivos sejam não-estacionários. 
# Para evitar esse problema, use a opção "method = c('ML')" dentro da função "Arima()" 
# Se for o caso (que não é o desse exemplo), pode ser também inserido o parâmetro "lambda = lambda" para receber o resultado "lambda = BoxCox.lambda(y)".

# Residuos
par(mfrow=c(1,1))
residuos <- fit2$residuals %>% window(start=c(1960,1))
par(mfrow=c(1,3))
plot(residuos,main="Resíduos após ininicialização do modelo")
qqnorm(residuos);qqline(residuos)
acf(residuos,lag.max=12*5)
shapiro <- shapiro.test(residuos)
kpss <- kpss.test(residuos)
box <- Box.test(residuos,lag=15,type='Ljung-Box')

# resumo dos testes
df.residuos <- data.frame("Teste"=c("Shapiro-Wilk","KPSS","Ljung-Box"),
                          "P-valor"=c(shapiro$p.value,kpss$p.value,box$p.value),
                          check.names=FALSE)
