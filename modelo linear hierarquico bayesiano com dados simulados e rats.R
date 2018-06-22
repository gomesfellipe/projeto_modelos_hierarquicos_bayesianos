# ############################################################################################################
# ############################################# Dados simulados ##############################################
# ############################################################################################################

# carregando dependencias

source("dependencies.R",encoding = "UTF-8")

# Amostrador de Gibbs para modelo linear bayesiano --------------------------

############################
#     Gerando a amostra
############################

set.seed(22-03-2018)

x             = c(8, 15, 22, 29, 36)
xbarra        = mean(x)
x             = x - xbarra
n             = 30
t             = 5
alpha_c       = 20
beta_c        = 2
alpha         = beta = NULL
taualpha      = 1/0.2
taubeta       = 1/0.2
y             = matrix(NA,n,t)
tau           = 1
e             = matrix(rnorm(n*t,0,sqrt(1/tau)),n,t)
for(i in 1:n){
  alpha[i]    = rnorm(1,alpha_c,sqrt(1/taualpha))
  beta[i]     = rnorm(1,beta_c,sqrt(1/taubeta))
  for(j in 1:t){
    y[i,j]    = alpha[i] + beta[i]*x[j] + e[i,j]
  }}

hist(y)

############################
#   Parametros a  Priori
############################

m.alpha = 0
V.alpha =1/0.0001

m.beta  =0
V.beta  =1/0.0001

a.tau   =0.001
b.tau   =0.001

a.alpha =0.001
b.alpha =0.001

a.beta  =0.001
b.beta  =0.001

############################
#    Valores da cadeia
############################

nsim              = 150000
burnin            =  nsim/2

#Criando os objetos:
cadeia.alpha.i    = matrix(NA,nsim,n)
cadeia.beta.i     = matrix(NA,nsim,n)
cadeia.tau.c      = matrix(NA,nsim,1)
cadeia.alpha.c    = matrix(NA,nsim,1)
cadeia.beta.c     = matrix(NA,nsim,1)
cadeia.tau.alpha  = matrix(NA,nsim,1)
cadeia.tau.beta   = matrix(NA,nsim,1)

#Chutes iniciais:
cadeia.alpha.i[1,]      = 0
cadeia.beta.i[1,]       = 0
cadeia.tau.c[1,]        = 1
cadeia.alpha.c[1,]      = 0
cadeia.beta.c[1,]       = 0
cadeia.tau.alpha[1,]    = 1
cadeia.tau.beta[1,]     = 1

# Obs.: Indice da cadeia sera k
# ----------------------------------------------------------------------------------

############################
#    Algoritimo da cadeia
############################

# Sementes:
set.seed(1)


#Criar uma barra de processo e acompanhar o carregamento:
pb <- txtProgressBar(min = 0, max = nsim, style = 3);antes = Sys.time()

for(k in 2:nsim){
  
  soma = 0
  #Cadeia alpha.i e beta.i
  for(i in 1:n){
    dccp.tau.alpha      = 1  /  ((t*cadeia.tau.c[k-1])  +  cadeia.tau.alpha[k-1])
    dccp.alpha.c        = (  cadeia.tau.c[k-1]*sum( y[i,] - cadeia.beta.i[k-1,i]*x  )  +  
                               cadeia.tau.alpha[k-1]*cadeia.alpha.c[k-1] )  * dccp.tau.alpha
    cadeia.alpha.i[k,i] = rnorm(1,dccp.alpha.c,sqrt(dccp.tau.alpha))
    
    
    dccp.tau.beta       = 1  /  ((cadeia.tau.c[k-1]*sum(x^2))  +  cadeia.tau.beta[k-1])
    u                   = (y[i,]-cadeia.alpha.i[k,i])*x
    dccp.beta.c         = (cadeia.tau.c[k-1]*sum(u) +  (cadeia.tau.beta[k-1]*cadeia.beta.c[k-1]))  *  dccp.tau.beta
    cadeia.beta.i[k,i]  = rnorm(1,dccp.beta.c,sqrt(dccp.tau.beta))
    
    for(j in 1:t){soma = soma + (y[i,j] - cadeia.alpha.i[k,i]-cadeia.beta.i[k,i]*x[j])^2}
  }
  
  # cadeia de tau
  dccp.a                 = ((n*t)/2) + a.tau 
  dccp.beta              = b.tau + 0.5  * soma
  cadeia.tau.c[k]        = rgamma(1,dccp.a,dccp.beta)
  #cadeia.tau.c[k]       = tau
  
  #cadeia de tau.alpha
  dccp.a.alpha           = (n/2) + a.alpha
  dccp.b.alpha           = b.alpha + 0.5 * sum((cadeia.alpha.i[k,]-cadeia.alpha.c[k-1,])^2) 
  cadeia.tau.alpha[k]    = rgamma(1,dccp.a.alpha,dccp.b.alpha)
  #cadeia.tau.alpha[k]   = taualpha
  
  #cadeia de tau.beta
  dccp.a.beta            = (n/2) + a.beta
  dccp.b.beta            = b.beta + 0.5 * sum((cadeia.beta.i[k,]-cadeia.beta.c[k-1,])^2) 
  cadeia.tau.beta[k]     = rgamma(1,dccp.a.beta,dccp.b.beta)
  #cadeia.tau.beta[k]    = taubeta
  
  #Cadeia de alpha.c
  dccp.V.alpha           = 1  /  (n*cadeia.tau.alpha[k]  +  1/V.alpha)
  dccp.m.alpha           = (cadeia.tau.alpha[k]  *  sum(cadeia.alpha.i[k,])  +  m.alpha/V.alpha) *  dccp.V.alpha
  cadeia.alpha.c[k]      = rnorm(1,dccp.m.alpha,sqrt(dccp.V.alpha))
  #cadeia.alpha.c[k]     = alpha_c
  
  #Cadeia de beta.c
  dccp.V.beta            = 1  /  (n*cadeia.tau.beta[k]  +  1/V.beta)
  dccp.m.beta            = (cadeia.tau.beta[k]  *  sum(cadeia.beta.i[k,])  +  m.beta/V.beta) *  dccp.V.beta
  cadeia.beta.c[k]       = rnorm(1,dccp.m.beta,sqrt(dccp.V.beta))
  #cadeia.beta.c[k]      = beta_c
  
  # update barra de processo
  setTxtProgressBar(pb, k)
  
}; close(pb); depois = Sys.time() - antes; depois #Encerrando barra de processo e tempo da cadeia


############################
#    Resultados da cadeia
############################

# Juntando resultados:
inds = seq(burnin,nsim,50)
results <- cbind(cadeia.alpha.c,
                 cadeia.beta.c,
                 cadeia.tau.c,
                 cadeia.tau.alpha,
                 cadeia.tau.beta) %>% as.data.frame() %>% .[inds,]

real    <- c(alpha_c,beta_c,tau,taualpha,taubeta)
name    <- c(expression(alpha[c]), expression(beta[c]), expression(tau[c]), expression(tau[alpha]), expression(tau[beta]))

# Cadeia
png("imagens/cadeia_hierarquica_sim.png",1100,1000)
cadeia(results,name,real)
dev.off()

# Histograma e densidade
png("imagens/hist_hierarquica_sim.png",830,610)
g1 <- hist_den(results[,1],name = name[1], p = real[1])
g2 <- hist_den(results[,2],name = name[2], p = real[2])
g3 <- hist_den(results[,3],name = name[3], p = real[3])
g4 <- hist_den(results[,4],name = name[4], p = real[4])
g5 <- hist_den(results[,5],name = name[5], p = real[5])
grid.arrange(g1,g2,g3,g4,g5,ncol=1)
dev.off()

# ACF
png("imagens/acf_hierarquica_sim.png",1000,610)
FAC(results)
dev.off()

# Coeficientes:
coeficientes <- coeficientes_hierarquico(cadeia.alpha.i[inds,],alpha,
                         cadeia.beta.i[inds,], beta,
                         cadeia.alpha.c, alpha_c,
                         cadeia.beta.c, beta_c,
                         cadeia.tau.c,tau,
                         cadeia.tau.alpha, taualpha,
                         cadeia.tau.beta, taubeta,
                         ncol(results)+60
)

# sumario:
sumario <- coeficientes[[1]];sumario
# credibilidade:
tab <- coeficientes[[2]];tab



tab %>% tail(5)  %>% .[,c(1,3,4)] %>% cbind(sd = results %>% apply(2,sd) %>% round(4)) %>%  xtable()

# Com GGPLOT para o texto -------------------------------------------------

##############
# Para alpha #
##############

png("imagens/caterplot_a_hierarquica_sim.png",800,600)
cbind(alpha=1:30,tab[1:30,])%>%
  as.data.frame%>%
  ggplot( aes(x=alpha, y=Média)) +  #, colour=supp, group=supp  <- separar por cadeia?
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), colour="black", width=.1, position=position_dodge(0.1)) +
  # geom_line(position=pd) +
  geom_point(aes(x=alpha,y=`Param. Pop.`),position=position_dodge(0.1), size=2, shape=21, fill="blue") + # 21 is filled circle
  geom_point(position=position_dodge(0.1), size=3, shape=21, fill="white") + # 21 is filled circle
  xlab("Cadeia de alpha") +
  ylab("Intervalo de credibilidade") +
  # scale_colour_hue(name="Supplement type",    # Legend label, use darker colors
  #                  breaks=c("OJ", "VC"),
  #                  labels=c("Orange juice", "Ascorbic acid"),
  #                  l=40) +                    # Use darker colors, lightness=40
  ggtitle(" ") +
  expand_limits(y=0) +                        # Expand y range
  scale_y_continuous(limits=c(min(tab[1:30,1]),max(tab[1:30,4])),breaks=seq(min(tab[1:30,1]),max(tab[1:30,4]),0.2)) +         # Set tick every 4
  theme_bw()+ 
  geom_hline(yintercept = alpha_c, linetype="dashed", color = "red",size=1.3)+theme(axis.text=element_text(size=12),
                                                                                    axis.title=element_text(size=14,face="bold"))
# +  # Para incluir legenda
#   theme(legend.justification=c(1,0),
#         legend.position=c(1,0))               # Position legend in bottom right
dev.off()

##############
# Para beta  #
##############

png("imagens/caterplot_b_hierarquica_sim.png",800,600)
cbind(beta=31:60,tab[31:60,])%>%
  as.data.frame%>%
  ggplot( aes(x=beta, y=Média)) +  #, colour=supp, group=supp  <- separar por cadeia?
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), colour="black", width=.1, position=position_dodge(0.1)) +
  geom_point(position=position_dodge(0.1), size=3, shape=21, fill="white") + # 21 is filled circle
  geom_point(aes(x=beta,y=`Param. Pop.`),position=position_dodge(0.1), size=2, shape=21, fill="blue") + # 21 is filled circle
  # geom_line(position=pd) +
  xlab("Cadeia de beta") +
  ylab("Intervalo de credibilidade") +
  # scale_colour_hue(name="Supplement type",    # Legend label, use darker colors
  #                  breaks=c("OJ", "VC"),
  #                  labels=c("Orange juice", "Ascorbic acid"),
  #                  l=40) +                    # Use darker colors, lightness=40
  ggtitle(" ") +
  expand_limits(y=0) +                        # Expand y range
  scale_y_continuous(limits=c(min(tab[31:60,1]),max(tab[31:60,4])),breaks=seq(min(tab[31:60,1]),max(tab[31:60,4]),0.2)) +         # Set tick every 4
  theme_bw()+ 
  geom_hline(yintercept = beta_c, linetype="dashed", color = "red",size=1.3)+theme(axis.text=element_text(size=12),
                                                                                   axis.title=element_text(size=14,face="bold"))
# +  # Para incluir legenda
#   theme(legend.justification=c(1,0),
#         legend.position=c(1,0))               # Position legend in bottom right
dev.off()

# ############################################################################################################
# #############################################   Dados reais   ##############################################
# ############################################################################################################

# Amostrador de Gibbs para modelo linear bayesiano --------------------------

############################
#     Amostra utilizada
############################

rats = read.table("rats.txt",header=T)

y=rats 
x=c(8, 15, 22, 29, 36)
xbarra=mean(x)
x=x-xbarra
n=nrow(y)
t=ncol(y)

############################
#   Parametros a  Priori
############################

m.alpha = 0
V.alpha =1/0.0001

m.beta  =0
V.beta  =1/0.0001

a.tau   =0.001
b.tau   =0.001

a.alpha =0.001
b.alpha =0.001

a.beta  =0.001
b.beta  =0.001

############################
#    Valores da cadeia
############################

nsim              = 50000
burnin            =  nsim/3

#Criando os objetos:
cadeia.alpha.i    = matrix(NA,nsim,n)
cadeia.beta.i     = matrix(NA,nsim,n)
cadeia.tau.c      = matrix(NA,nsim,1)
cadeia.alpha.c    = matrix(NA,nsim,1)
cadeia.beta.c     = matrix(NA,nsim,1)
cadeia.tau.alpha  = matrix(NA,nsim,1)
cadeia.tau.beta   = matrix(NA,nsim,1)

#Chutes iniciais:
cadeia.alpha.i[1,]      = 0
cadeia.beta.i[1,]       = 0
cadeia.tau.c[1,]        = 1
cadeia.alpha.c[1,]      = 0
cadeia.beta.c[1,]       = 0
cadeia.tau.alpha[1,]    = 1
cadeia.tau.beta[1,]     = 1

# Obs.: Indice da cadeia sera k
# ----------------------------------------------------------------------------------

############################
#    Algoritimo da cadeia
############################

# Sementes:
set.seed(1)


#Criar uma barra de processo e acompanhar o carregamento:
pb <- txtProgressBar(min = 0, max = nsim, style = 3);antes = Sys.time()

for(k in 2:nsim){
  
  soma = 0
  #Cadeia alpha.i e beta.i
  for(i in 1:n){
    dccp.tau.alpha      = 1  /  ((t*cadeia.tau.c[k-1])  +  cadeia.tau.alpha[k-1])
    dccp.alpha.c        = (  cadeia.tau.c[k-1]*sum( y[i,] - cadeia.beta.i[k-1,i]*x  )  +  
                               cadeia.tau.alpha[k-1]*cadeia.alpha.c[k-1] )  * dccp.tau.alpha
    cadeia.alpha.i[k,i] = rnorm(1,dccp.alpha.c,sqrt(dccp.tau.alpha))
    
    
    dccp.tau.beta       = 1  /  ((cadeia.tau.c[k-1]*sum(x^2))  +  cadeia.tau.beta[k-1])
    u                   = (y[i,]-cadeia.alpha.i[k,i])*x
    dccp.beta.c         = (cadeia.tau.c[k-1]*sum(u) +  (cadeia.tau.beta[k-1]*cadeia.beta.c[k-1]))  *  dccp.tau.beta
    cadeia.beta.i[k,i]  = rnorm(1,dccp.beta.c,sqrt(dccp.tau.beta))
    
    for(j in 1:t){soma = soma + (y[i,j] - cadeia.alpha.i[k,i]-cadeia.beta.i[k,i]*x[j])^2}
  }
  
  # cadeia de tau
  dccp.a                 = ((n*t)/2) + a.tau 
  dccp.beta              = b.tau + 0.5  * soma
  cadeia.tau.c[k]        = rgamma(1,dccp.a,dccp.beta)
  #cadeia.tau.c[k]       = tau
  
  #cadeia de tau.alpha
  dccp.a.alpha           = (n/2) + a.alpha
  dccp.b.alpha           = b.alpha + 0.5 * sum((cadeia.alpha.i[k,]-cadeia.alpha.c[k-1,])^2) 
  cadeia.tau.alpha[k]    = rgamma(1,dccp.a.alpha,dccp.b.alpha)
  #cadeia.tau.alpha[k]   = taualpha
  
  #cadeia de tau.beta
  dccp.a.beta            = (n/2) + a.beta
  dccp.b.beta            = b.beta + 0.5 * sum((cadeia.beta.i[k,]-cadeia.beta.c[k-1,])^2) 
  cadeia.tau.beta[k]     = rgamma(1,dccp.a.beta,dccp.b.beta)
  #cadeia.tau.beta[k]    = taubeta
  
  #Cadeia de alpha.c
  dccp.V.alpha           = 1  /  (n*cadeia.tau.alpha[k]  +  1/V.alpha)
  dccp.m.alpha           = (cadeia.tau.alpha[k]  *  sum(cadeia.alpha.i[k,])  +  m.alpha/V.alpha) *  dccp.V.alpha
  cadeia.alpha.c[k]      = rnorm(1,dccp.m.alpha,sqrt(dccp.V.alpha))
  #cadeia.alpha.c[k]     = alpha_c
  
  #Cadeia de beta.c
  dccp.V.beta            = 1  /  (n*cadeia.tau.beta[k]  +  1/V.beta)
  dccp.m.beta            = (cadeia.tau.beta[k]  *  sum(cadeia.beta.i[k,])  +  m.beta/V.beta) *  dccp.V.beta
  cadeia.beta.c[k]       = rnorm(1,dccp.m.beta,sqrt(dccp.V.beta))
  #cadeia.beta.c[k]      = beta_c
  
  # update barra de processo
  setTxtProgressBar(pb, k)
  
}; close(pb); depois = Sys.time() - antes; depois #Encerrando barra de processo e tempo da cadeia


############################
#    Resultados da cadeia
############################

# Juntando resultados:
inds = seq(burnin,nsim,50)
results <- cbind(cadeia.alpha.c,
                 cadeia.beta.c,
                 cadeia.tau.c,
                 cadeia.tau.alpha,
                 cadeia.tau.beta) %>% as.data.frame() %>% .[inds,]

name    <- c(expression(alpha[c]), expression(beta[c]), expression(tau[c]), expression(tau[alpha]), expression(tau[beta]))

# Cadeia
png("imagens/cadeia_hierarquica_rats.png",830,610)
cadeia(results,name)
dev.off()

# Histograma e densidade
png("imagens/hist_hierarquica_rats.png",830,610)
g1 <- hist_den(results[,1],name = name[1])
g2 <- hist_den(results[,2],name = name[2])
g3 <- hist_den(results[,3],name = name[3])
g4 <- hist_den(results[,4],name = name[4])
g5 <- hist_den(results[,5],name = name[5])
grid.arrange(g1,g2,g3,g4,g5,ncol=1)
dev.off()

# ACF
png("imagens/acf_hierarquica_rats.png",1000,610)
FAC(results)
dev.off()

###################
# continuar aqui
##################

# Coeficientes:
coeficientes <- coeficientes_hierarquico2(cadeia.alpha.i[inds,],
                                         cadeia.beta.i[inds,],
                                         cadeia.alpha.c,
                                         cadeia.beta.c, 
                                         cadeia.tau.c,
                                         cadeia.tau.alpha,
                                         cadeia.tau.beta,
                                         ncol(results)+60
)

# sumario:
sumario <- coeficientes[[1]];sumario
# credibilidade:
tab <- coeficientes[[2]];tab



tab %>% tail(5)  %>% cbind(sd = results %>% apply(2,sd) %>% round(4)) %>%  xtable()

# Com GGPLOT para o texto -------------------------------------------------

##############
# Para alpha #
##############

png("imagens/caterplot_a_hierarquica_rats.png",800,600)
cbind(alpha=1:30,tab[1:30,])%>%
  as.data.frame%>%
  ggplot( aes(x=alpha, y=Média)) +  #, colour=supp, group=supp  <- separar por cadeia?
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), colour="black", width=.1, position=position_dodge(0.1)) +
  # geom_line(position=pd) +
  geom_point(aes(x=alpha,y=Média),position=position_dodge(0.1), size=2, shape=21, fill="blue") + # 21 is filled circle
  geom_point(position=position_dodge(0.1), size=3, shape=21, fill="white") + # 21 is filled circle
  xlab("Cadeia de alpha") +
  ylab("Intervalo de credibilidade") +
  # scale_colour_hue(name="Supplement type",    # Legend label, use darker colors
  #                  breaks=c("OJ", "VC"),
  #                  labels=c("Orange juice", "Ascorbic acid"),
  #                  l=40) +                    # Use darker colors, lightness=40
  ggtitle(" ") +
  expand_limits(y=0) +                        # Expand y range
  scale_y_continuous(limits=c(min(tab[1:30,1]),max(tab[1:30,3]))) +         # Set tick every 4
  scale_x_continuous(breaks = seq(1,30,1)) +         # Set tick every 4
  theme_bw()
# + 
  # geom_hline(yintercept = alpha_c, linetype="dashed", color = "red",size=1.3)+theme(axis.text=element_text(size=12),
                                                                                    # axis.title=element_text(size=14,face="bold"))
# +  # Para incluir legenda
#   theme(legend.justification=c(1,0),
#         legend.position=c(1,0))               # Position legend in bottom right
dev.off()

##############
# Para beta  #
##############

png("imagens/caterplot_b_hierarquica_rats.png",800,600)
cbind(beta=1:30,tab[31:60,])%>%
  as.data.frame%>%
  ggplot( aes(x=beta, y=Média)) +  #, colour=supp, group=supp  <- separar por cadeia?
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), colour="black", width=.1, position=position_dodge(0.1)) +
  # geom_line(position=pd) +
  geom_point(aes(x=beta,y=Média),position=position_dodge(0.1), size=2, shape=21, fill="blue") + # 21 is filled circle
  geom_point(position=position_dodge(0.1), size=3, shape=21, fill="white") + # 21 is filled circle
  xlab("Cadeia de beta") +
  ylab("Intervalo de credibilidade") +
  # scale_colour_hue(name="Supplement type",    # Legend label, use darker colors
  #                  breaks=c("OJ", "VC"),
  #                  labels=c("Orange juice", "Ascorbic acid"),
  #                  l=40) +                    # Use darker colors, lightness=40
  ggtitle(" ") +
  expand_limits(y=0) +                        # Expand y range
  scale_y_continuous(limits=c(min(tab[31:60,1]),max(tab[31:60,3]))) +         # Set tick every 4
  scale_x_continuous(breaks = seq(1,30,1)) +         # Set tick every 4
  theme_bw()
# + 
# geom_hline(yintercept = beta_c, linetype="dashed", color = "red",size=1.3)+theme(axis.text=element_text(size=12),
# axis.title=element_text(size=14,face="bold"))
# +  # Para incluir legenda
#   theme(legend.justification=c(1,0),
#         legend.position=c(1,0))               # Position legend in bottom right
dev.off()
tab
