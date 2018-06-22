# ############################################################################################################
# ############################################# Dados simulados ##############################################
# ############################################################################################################

# carregando dependencias

source("dependencies.R",encoding = "UTF-8")

############################
#     Gerando a amostra
############################

set.seed(12)

# Amostra que sera utilizada:
# N=1000, \beta_0 = 1, \beta_1 = 0,5, \tau = 2 e X_i ~ N(0,1)
n   <- 1000
b0  <- 1
b1  <- 0.5
tau <- 2
x   <- rnorm(n)
y   <- b0 + b1 * x + rnorm(n,0,sqrt(1/tau))


#Conferindo normalidade dos dados
shapiro.test(y)
shapiro.test(x)


############################
#   Parametros a  Priori
############################

# \bs{\theta} = (\beta_0, \beta_1, \tau)
# m_0 = m_1 = 0, \sigma_0^2 = sigma_1^2 = 100, a=0,1 e b=0,1
#Parametros para b0 ~ N(mu0, sig0)
mu0   = 0
sig0  = 1000


#Parametros para b1 ~ N(mu1, sig1)
mu1   = 0
sig1  = 1000

#Parametros para tau ~ G(a,b)
a     = 0.1
b     = 0.1


############################
#    Valores da cadeia
############################

nsim          = 3*10000
cadeia.b0     = rep(0,nsim)
cadeia.b1     = rep(0,nsim)
cadeia.tau    = rep(0,nsim)

# #Chutes iniciais: 

cadeia.b0[1]   = 0
cadeia.b1[1]   = 0
cadeia.tau[1]  = 1

############################
#    Algoritimo da cadeia
############################

pb <- txtProgressBar(min = 0, max = nsim, style = 3) # iniciando barra de processo
for (k in 2:nsim){
  
  #Cadeia tau
  cadeia.tau[k]     = rgamma(1, (n/2)+a, b + (sum(( y- cadeia.b0[k-1] - (cadeia.b1[k-1]*x)  )^2)/2) )
  
  # Cadeia B0
  c0                = (n*cadeia.tau[k]) + (1/sig0)
  m0                = ( cadeia.tau[k]*sum(y) - (cadeia.tau[k]*cadeia.b1[k-1]*sum(x)) + (mu0/sig0)  )   / c0
  cadeia.b0[k]      = rnorm(1, m0, 1/sqrt(c0))
  
  # Cadeia B1
  c1                =  ( sum(x^2)*cadeia.tau[k] ) + (1/sig1)
  m1                =  ( (cadeia.tau[k]*sum(x*y)) - (cadeia.tau[k]*cadeia.b0[k]*sum(x))  + (mu1/sig1)  )   /  c1
  cadeia.b1[k]      = rnorm(1, m1, 1/sqrt(c1))
  
  
  setTxtProgressBar(pb, k)
  
}; close(pb) #Encerrando barra de processo

############################
#    Resultados da cadeia
############################

# Juntando resultados:
inds    <- seq(nsim/2,nsim) # Definindo os indices
results <- cbind(cadeia.b0,cadeia.b1,cadeia.tau) %>% as.data.frame() %>% .[inds,]
real    <- c(b0,b1,tau)
name    <- c(expression(beta[0]),expression(beta[1]),expression(tau))

# Cadeia
cadeia(results,name,real)

# Histograma e densidade
g1 <- hist_den(results[,1],name = name[1], p = real[1])
g2 <- hist_den(results[,2],name = name[2], p = real[2])
g3 <- hist_den(results[,3],name = name[3], p = real[3])
grid.arrange(g1,g2,g3,ncol=1)

# ACF
FAC(results)

# Resultados:
coeficientes(results)              # Sem parametros reais
coeficientes(results,real = real)  # Com parametro real


# Reta do modelo classico -------------------------------------------------

plot(x,y)
modelo.classico=lm(y~ 1 + x)
a.classico     =modelo.classico$coefficients[1]
b.classico     =modelo.classico$coefficients[2]
abline(a       =a.classico,b=b.classico,col="blue")

# Reta do modelo bayesiano
plot(x,y)
a.bayes =mean(results[,1])
b.bayes =mean(results[,2])
abline(a=a.bayes,b=b.bayes,col="red")


# Com GGPLOT para o texto -------------------------------------------------

# pacotes=c("stringr", "ggplot2", "ggExtra")
# install.packages(pacotes)

library(stringr)
library(ggplot2)
library(ggExtra)

# Texto da imagem
text.classico=str_c("Modelo Classico: ","y = ",round(a.classico,4)," x + ",round(b.classico,4))
text.bayes=str_c("Modelo Bayesiano: ","y = ",round(a.bayes,4)," x + ",round(b.bayes,4))

#Gerando o scatter.plot
cbind(y,x) %>% 
  as.data.frame %>% 
  ggplot(aes(y=y, x=x))+
  geom_point()+
  geom_smooth(method = "lm",se=F,col="red")+
  theme_classic()+
  geom_abline( slope = b.bayes,intercept = a.bayes,col="blue")+
  labs(title="",
       x="Covariável",
       y="Reposta")
  # scale_y_continuous(breaks = seq(0,30,5),limits = c(0,30))+
  # scale_x_continuous(breaks = seq(0,130,10),limits = c(0,130))+
  # annotate("text", label=text.classico,  x=100, y=5)+
  # annotate("text", label=text.bayes,  x=101, y=3)

# ############################################################################################################
# #############################################   Dados reais   ##############################################
# ############################################################################################################

# Amostrador de Gibbs para modelo linear bayesiano --------------------------

############################
#     Amostra utilizada
############################

y   = cars$speed
x   = cars$dist
n   = nrow(cars)

#Conferindo normalidade dos dados
shapiro.test(y)
shapiro.test(x)

############################
#   Parametros a  Priori
############################

# \bs{\theta} = (\beta_0, \beta_1, \tau)
# m_0 = m_1 = 0, \sigma_0^2 = sigma_1^2 = 100, a=0,1 e b=0,1
#Parametros para b0 ~ N(mu0, sig0)
mu0   = 0
sig0  = 100

#Parametros para b1 ~ N(mu1, sig1)
mu1   = 0
sig1  = 1000

#Parametros para tau ~ G(a,b)
a     = 2
b     = 2

############################
#    Valores da cadeia
############################

nsim          = 3*10000
cadeia.b0     = rep(0,nsim)
cadeia.b1     = rep(0,nsim)
cadeia.tau    = rep(0,nsim)



# #Chutes iniciais:

cadeia.b0[1]   = 0
cadeia.b1[1]   = 0
cadeia.tau[1]  = 1

############################
#    Algoritimo da cadeia
############################

pb <- txtProgressBar(min = 0, max = nsim, style = 3)
for (k in 2:nsim){
  
  #Cadeia tau
  cadeia.tau[k]     = rgamma(1, (n/2)+a, b + (sum(( y- cadeia.b0[k-1] - (cadeia.b1[k-1]*x)  )^2)/2) )
  
  # Cadeia B0
  c0                = (n*cadeia.tau[k]) + (1/sig0)
  m0                = ( cadeia.tau[k]*sum(y) - (cadeia.tau[k]*cadeia.b1[k-1]*sum(x)) + (mu0/sig0)  )   / c0
  cadeia.b0[k]      = rnorm(1, m0, 1/sqrt(c0))
  
  # Cadeia B1
  c1                =  ( sum(x^2)*cadeia.tau[k] ) + (1/sig1)
  m1                =  ( (cadeia.tau[k]*sum(x*y)) - (cadeia.tau[k]*cadeia.b0[k]*sum(x))  + (mu1/sig1)  )   /  c1
  cadeia.b1[k]      = rnorm(1, m1, 1/sqrt(c1))
  
  
  setTxtProgressBar(pb, k)
  
}; close(pb) #Encerrando barra de processo


############################
#    Resultados da cadeia
############################

# Juntando resultados:
inds     <- seq(nsim/2,nsim) # Definindo os indices
results  <- cbind(cadeia.b0,cadeia.b1,cadeia.tau) %>% as.data.frame() %>% .[inds,]
classico <- c(coefficients(lm(cars)),1/var(lm(cars)$residuals))
name     <- c(expression(beta[0]),expression(beta[1]),expression(tau))

# Cadeia
cadeia(results,name)

# Histograma e densidade
g1 <- hist_den(results[,1],name = name[1])
g2 <- hist_den(results[,2],name = name[2])
g3 <- hist_den(results[,3],name = name[3])
grid.arrange(g1,g2,g3,ncol=1)

# ACF
FAC(results)

# Resultados:
coeficientes(results)              # Sem parametros reais
coeficientes(results,real = classico)  # Com parametro classico


# Reta do modelo classico -------------------------------------------------

plot(x,y)
modelo.classico=lm(y~ 1 + x)
a.classico     =modelo.classico$coefficients[1]
b.classico     =modelo.classico$coefficients[2]
abline(a       =a.classico,b=b.classico,col="blue")

# Reta do modelo bayesiano
plot(x,y)
a.bayes =mean(results[,1])
b.bayes =mean(results[,2])
abline(a=a.bayes,b=b.bayes,col="red")


# Com GGPLOT para o texto -------------------------------------------------

# Texto da imagem
text.classico=str_c("Modelo Classico: ","y = ",round(a.classico,4)," x + ",round(b.classico,4))
text.bayes   =str_c("Modelo Bayesiano: ","y = ",round(a.bayes,4)," x + ",round(b.bayes,4))

#Gerando o scatter.plot
cbind(y, x) %>%
  as.data.frame %>%
  ggplot(aes(y = y, x = x)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, col = "red") +
  theme_classic() +
  geom_abline(slope = b.bayes,
              intercept = a.bayes,
              col = "blue") +
  labs(title = "Relação entre a Distância e a Velocidade com \nreta do modelo linear clássico vs bayesiano",
       x = "Distância",
       y = "Velocidade")
# scale_y_continuous(breaks = seq(0,30,5),limits = c(0,30))+
# scale_x_continuous(breaks = seq(0,130,10),limits = c(0,130))+
# annotate("text", label=text.classico,  x=100, y=5)+
# annotate("text", label=text.bayes,  x=101, y=3)
