
# Pacotes -----------------------------------------------------------------

packages = c('dplyr', 'ggplot2','gridExtra', 'dplyr','purrr', 'xtable','stringr', 'ggExtra','plotly','magrittr')

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    suppressMessages(library(package, character.only=T))
  }
}


# Amostrador de Gibbs -------------------------------------------------------

# As 3 full conditionals fecham por conjugacao: tau | resto ~ Gama (conjugada
# com a verossimilhanca Normal) e b0, b1 | resto ~ Normal (conjugada Normal-
# Normal). Cada iteracao amostra as 3 nessa ordem, condicionando sempre no
# valor mais recente dos outros dois (Gibbs).
gibbs_regressao_linear <- function(y, x, n, nsim, mu0, sig0, mu1, sig1, a, b){

  cadeia.b0     = rep(0,nsim)
  cadeia.b1     = rep(0,nsim)
  cadeia.tau    = rep(0,nsim)

  # Chutes iniciais:
  cadeia.b0[1]   = 0
  cadeia.b1[1]   = 0
  cadeia.tau[1]  = 1

  pb <- txtProgressBar(min = 0, max = nsim, style = 3) # iniciando barra de processo
  for (k in 2:nsim){

    # tau | b0, b1, y ~ Gama(n/2 + a, b + SQE/2)
    cadeia.tau[k]     = rgamma(1, (n/2)+a, b + (sum(( y- cadeia.b0[k-1] - (cadeia.b1[k-1]*x)  )^2)/2) )

    # b0 | b1, tau, y ~ Normal
    c0                = (n*cadeia.tau[k]) + (1/sig0)
    m0                = ( cadeia.tau[k]*sum(y) - (cadeia.tau[k]*cadeia.b1[k-1]*sum(x)) + (mu0/sig0)  )   / c0
    cadeia.b0[k]      = rnorm(1, m0, 1/sqrt(c0))

    # b1 | b0, tau, y ~ Normal
    c1                =  ( sum(x^2)*cadeia.tau[k] ) + (1/sig1)
    m1                =  ( (cadeia.tau[k]*sum(x*y)) - (cadeia.tau[k]*cadeia.b0[k]*sum(x))  + (mu1/sig1)  )   /  c1
    cadeia.b1[k]      = rnorm(1, m1, 1/sqrt(c1))

    setTxtProgressBar(pb, k)

  }; close(pb) #Encerrando barra de processo

  cbind(cadeia.b0,cadeia.b1,cadeia.tau) %>% as.data.frame()
}

# Full conditionals fecham por conjugacao em 2 niveis: alpha_i, beta_i | resto
# ~ Normal (nivel do individuo i), tau, tau.alpha, tau.beta | resto ~ Gama
# (precisoes), alpha_c, beta_c | resto ~ Normal (nivel da populacao). Prefixo
# dccp. = "distribuicao condicional completa a posteriori", o parametro da
# full conditional amostrada em seguida.
gibbs_hierarquico <- function(y, x, n, t, nsim,
                               m.alpha, V.alpha,
                               m.beta, V.beta,
                               a.tau, b.tau,
                               a.alpha, b.alpha,
                               a.beta, b.beta){

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

  #Criar uma barra de processo e acompanhar o carregamento:
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)

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

    # tau | resto ~ Gama, usa soma dos quadrados dos residuos de todo i,j
    dccp.a                 = ((n*t)/2) + a.tau
    dccp.beta              = b.tau + 0.5  * soma
    cadeia.tau.c[k]        = rgamma(1,dccp.a,dccp.beta)

    # tau.alpha | resto ~ Gama (precisao dos alpha_i em torno de alpha_c)
    dccp.a.alpha           = (n/2) + a.alpha
    dccp.b.alpha           = b.alpha + 0.5 * sum((cadeia.alpha.i[k,]-cadeia.alpha.c[k-1,])^2)
    cadeia.tau.alpha[k]    = rgamma(1,dccp.a.alpha,dccp.b.alpha)

    # tau.beta | resto ~ Gama (precisao dos beta_i em torno de beta_c)
    dccp.a.beta            = (n/2) + a.beta
    dccp.b.beta            = b.beta + 0.5 * sum((cadeia.beta.i[k,]-cadeia.beta.c[k-1,])^2)
    cadeia.tau.beta[k]     = rgamma(1,dccp.a.beta,dccp.b.beta)

    # alpha_c | resto ~ Normal (media populacional dos alpha_i)
    dccp.V.alpha           = 1  /  (n*cadeia.tau.alpha[k]  +  1/V.alpha)
    dccp.m.alpha           = (cadeia.tau.alpha[k]  *  sum(cadeia.alpha.i[k,])  +  m.alpha/V.alpha) *  dccp.V.alpha
    cadeia.alpha.c[k]      = rnorm(1,dccp.m.alpha,sqrt(dccp.V.alpha))

    # beta_c | resto ~ Normal (media populacional dos beta_i)
    dccp.V.beta            = 1  /  (n*cadeia.tau.beta[k]  +  1/V.beta)
    dccp.m.beta            = (cadeia.tau.beta[k]  *  sum(cadeia.beta.i[k,])  +  m.beta/V.beta) *  dccp.V.beta
    cadeia.beta.c[k]       = rnorm(1,dccp.m.beta,sqrt(dccp.V.beta))

    # update barra de processo
    setTxtProgressBar(pb, k)

  }; close(pb) #Encerrando barra de processo

  list(alpha.i   = cadeia.alpha.i,
       beta.i    = cadeia.beta.i,
       tau.c     = cadeia.tau.c,
       alpha.c   = cadeia.alpha.c,
       beta.c    = cadeia.beta.c,
       tau.alpha = cadeia.tau.alpha,
       tau.beta  = cadeia.tau.beta)
}


# Cadeia ------------------------------------------------------------------

cadeia = function(df,name,p=NULL){
  if(is.null(p)){
    par(mar = c(5,5,2,2),
        mfrow = c(ncol(df),1))
    for(i in 1:ncol(df)){
      plot(df[,i],type="l",ylab=name[i],xlab="Iterações", cex.lab=1.8,cex.main=2)
      abline(h=quantile(df[,i],0.025),lty=3,col="blue",lwd=2)
      abline(h=quantile(df[,i],0.975),lty=3,col="blue",lwd=2)
    }
    par(mfrow=c(1,1))
  }else{
    par(mar = c(5,5,2,2),
        mfrow = c(ncol(df),1))
    for(i in 1:ncol(df)){
      plot(df[,i],type="l",ylab=name[i],xlab="Iterações", cex.lab=1.8,cex.main=2)
      abline(h=mean(p[i]), col="red",lwd=2)
      abline(h=quantile(df[,i],0.025),lty=3,col="blue",lwd=2)
      abline(h=quantile(df[,i],0.975),lty=3,col="blue",lwd=2)
    }
  }
}


# Histograma e densidade --------------------------------------------------

hist_den <- function(df, name = name,p=NULL){
  if(is.null(p)){
    g <- ggplot(data=df[1:length(df)/3] %>% as.data.frame(), aes(x=.)) + 
      geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                     colour="black", fill="white") +
      geom_density(alpha=.2, fill="lightgrey")+  # Overlay with transparent density plot
      labs(y="Densidade", x=name )+
      theme(axis.text.x = element_text(size=50))+
      theme(axis.text.y = element_text(size=50))+ 
      theme(axis.title.y = element_text(size = rel(1.8)))+
      theme_classic()+
      geom_density(data=df[(length(df)/3):(2*(length(df)/3))] %>% as.data.frame(),aes(x=., colour=I("darkred")))+
      geom_density(data=df[(2*(length(df)/3)):length(df)] %>% as.data.frame(),aes(x=., colour=I("darkblue")))
    
    return(g)
  }else{
    g <- ggplot(data=df[1:length(df)/3] %>% as.data.frame(), aes(x=.)) + 
      geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                     colour="black", fill="white") +
      geom_density(alpha=.2, fill="lightgrey")+  # Overlay with transparent density plot
      labs(y="Densidade", x=name )+
      theme(axis.text.x = element_text(size=50))+
      theme(axis.text.y = element_text(size=50))+ 
      theme(axis.title.y = element_text(size = rel(1.8)))+
      theme_classic()+
      geom_density(data=df[(length(df)/3):(2*(length(df)/3))] %>% as.data.frame(),aes(x=., colour=I("darkred")))+
      geom_density(data=df[(2*(length(df)/3)):length(df)] %>% as.data.frame(),aes(x=., colour=I("darkblue")))+ 
      geom_vline(xintercept = p,colour="red")
    
    return(g)
  }
  
}
  

# Autocorrelacao ----------------------------------------------------------

FAC <- function(df) {
  par(mar = c(5,5.5,2,2),
      mfrow = c(ncol(df),1))
  for(i in 1:ncol(df)){
    acf(x = df[,i],xlab=paste("Defasagem"),ylab="FAC",main="",
        cex.lab = 1.8,
        cex.main = 2
    )
  }
  par(mfrow=c(1,1))
}

#  Resultados dos coeficientes  -------------------------------------------

coeficientes <- function(df,real=NULL){
  if(is.null(real)){
  cbind(
  "Média" = apply(df,2,mean),
  "Desv. Pad." = apply(df,2,sd),
  "IC inf" = apply(df,2,function(x)quantile(x,0.025)),
  "IC sup" = apply(df,2,function(x)quantile(x,0.975))) %>% 
    round(4)
  }else{
    cbind(
      "Média" = apply(df,2,mean),
      "Desv. Pad." = apply(df,2,sd),
      "IC inf" = apply(df,2,function(x)quantile(x,0.025)),
      "IC sup" = apply(df,2,function(x)quantile(x,0.975)),
      "Real" = as.numeric(real)) %>% 
      round(4) 
    }
}

# Sumario do modelo hierárquico bayesiano ---------------------------------

coeficientes_hierarquico <- function(cadeia.alpha.i,alpha,
                                     cadeia.beta.i,beta,
                                     cadeia.alpha.c,alpha_c,
                                     cadeia.beta.c,beta_c,
                                     cadeia.tau.c,tau,
                                     cadeia.tau.alpha,taualpha,
                                     cadeia.tau.beta,taubeta,
                                     n){
  
  #Sumario do resultado dos parâmetros:
  sumario <- 
    cbind(apply(cadeia.alpha.i, 2, summary),
        apply(cadeia.beta.i, 2, summary)) %>%
        t %>%
        as.data.frame(row.names = c(paste(rep("alpha", n), 1:n), paste(rep("beta", n), 1:n))) %>% 
        rbind(
          cbind(
            cadeia.alpha.c,
            cadeia.beta.c,
            cadeia.tau.c,
            cadeia.tau.alpha,
            cadeia.tau.beta) %>% 
            as.data.frame %>% 
            apply(2,function(x) (summary(x) )) %>% 
            t %>% 
            as.data.frame(row.names = c("alpha.c", "beta.c", "tau.c", "tau.alpha", "tau.beta"))
        ) %>% 
        round(4)
  
      
#Quantos intervalos de credibilidade de 95% contem o real parametro: 
  intervalo.credibilidade <- 
    cbind(apply(cadeia.alpha.i, 2, function(x) quantile(x, probs = c(0.025, 0.975)) ),
          apply(cadeia.beta.i, 2, function(x) quantile(x, probs = c(0.025, 0.975)))) %>%
    t %>%
    as.data.frame(row.names = c(paste(rep("alpha", n), 1:n), paste(rep("beta", n), 1:n))) %>% 
    rbind(
      cbind(
        cadeia.alpha.c,
        cadeia.beta.c,
        cadeia.tau.c,
        cadeia.tau.alpha,
        cadeia.tau.beta) %>% 
        as.data.frame %>% 
        apply(2,function(x) quantile(x, probs = c(0.025, 0.975))) %>% 
        t %>% 
        as.data.frame(row.names = c("alpha.c", "beta.c", "tau.c", "tau.alpha", "tau.beta"))
    ) %>% 
    round(4)
    
    
  
  for(i in 1:length(alpha)){
    intervalo.credibilidade[i,3] = alpha[i]
    intervalo.credibilidade[i,4] = ifelse(alpha[i]>=intervalo.credibilidade[i,1] | alpha[i]<=intervalo.credibilidade[i,2],"Sim","Nao" )   
  }
  for(i in (length(alpha)+1):(length(alpha)+length(beta))){
    intervalo.credibilidade[i,3] = beta[i-length(alpha)]
    intervalo.credibilidade[i,4] = ifelse(beta[i-length(alpha)]>=intervalo.credibilidade[i,1] | beta[i-length(alpha)]<=intervalo.credibilidade[i,2],"Sim","Nao" )
  }
  intervalo.credibilidade["alpha.c",3:4]   = cbind(alpha_c,ifelse(alpha_c>=intervalo.credibilidade["alpha.c",1] | alpha_c<=intervalo.credibilidade["alpha.c",2],"Sim","Nao" ))
  intervalo.credibilidade["beta.c",3:4]    = cbind(beta_c,ifelse(beta_c>=intervalo.credibilidade["beta.c",1] | beta_c<=intervalo.credibilidade["beta.c",2],"Sim","Nao" ))
  intervalo.credibilidade["tau.c",3:4]     = cbind(tau,ifelse(tau>=intervalo.credibilidade["tau.c",1] | tau<=intervalo.credibilidade["tau.c",2],"Sim","Nao" ))
  intervalo.credibilidade["tau.alpha",3:4] = cbind(taualpha,ifelse(taualpha>=intervalo.credibilidade["tau.alpha",1] | taualpha<=intervalo.credibilidade["tau.alpha",2],"Sim","Nao" ))
  intervalo.credibilidade["tau.beta",3:4]  = cbind(taubeta,ifelse(taubeta>=intervalo.credibilidade["tau.beta",1] | taubeta<=intervalo.credibilidade["tau.beta",2],"Sim","Nao" ))
  
  # Manipulando dataframe para arredondar os resultados:
  intervalo.credibilidade[, -ncol(intervalo.credibilidade)] = apply(intervalo.credibilidade[, -ncol(intervalo.credibilidade)], 2, as.numeric)
  intervalo.credibilidade[, -ncol(intervalo.credibilidade)] = round(intervalo.credibilidade[, -ncol(intervalo.credibilidade)], 3)
  colnames(intervalo.credibilidade)[3:4] = c("Param. Pop.", "Estimou?")
  
  tab = as.data.frame(
    cbind(
      "2.5%" = intervalo.credibilidade[, 1],
      "Mediana" = sumario$Median,
      "Média" = sumario$Mean,
      "97.5%" = intervalo.credibilidade[, 2],
      intervalo.credibilidade[, 3:4]
    ), row.names = row.names(intervalo.credibilidade)
  )
  
  return(resultados = list(sumario, tab))
}

# coeficientes sem valor real ---------------------------------------------


coeficientes_hierarquico2 <- function(cadeia.alpha.i,
                                     cadeia.beta.i,
                                     cadeia.alpha.c,
                                     cadeia.beta.c,
                                     cadeia.tau.c,
                                     cadeia.tau.alpha,
                                     cadeia.tau.beta,
                                     n){
  
  #Sumario do resultado dos parâmetros:
  sumario <- 
    cbind(apply(cadeia.alpha.i, 2, summary),
        apply(cadeia.beta.i, 2, summary)) %>%
        t %>%
        as.data.frame(row.names = c(paste(rep("alpha", n), 1:n), paste(rep("beta", n), 1:n))) %>% 
        rbind(
          cbind(
            cadeia.alpha.c,
            cadeia.beta.c,
            cadeia.tau.c,
            cadeia.tau.alpha,
            cadeia.tau.beta) %>% 
            as.data.frame %>% 
            apply(2,function(x) (summary(x) )) %>% 
            t %>% 
            as.data.frame(row.names = c("alpha.c", "beta.c", "tau.c", "tau.alpha", "tau.beta"))
        ) %>% 
        round(4)
  
      
#Quantos intervalos de credibilidade de 95% contem o real parametro: 
  intervalo.credibilidade <- 
    cbind(apply(cadeia.alpha.i, 2, function(x) quantile(x, probs = c(0.025, 0.975)) ),
          apply(cadeia.beta.i, 2, function(x) quantile(x, probs = c(0.025, 0.975)))) %>%
    t %>%
    as.data.frame(row.names = c(paste(rep("alpha", n), 1:n), paste(rep("beta", n), 1:n))) %>% 
    rbind(
      cbind(
        cadeia.alpha.c,
        cadeia.beta.c,
        cadeia.tau.c,
        cadeia.tau.alpha,
        cadeia.tau.beta) %>% 
        as.data.frame %>% 
        apply(2,function(x) quantile(x, probs = c(0.025, 0.975))) %>% 
        t %>% 
        as.data.frame(row.names = c("alpha.c", "beta.c", "tau.c", "tau.alpha", "tau.beta"))
    ) %>% 
    round(4)
    
    
tab = as.data.frame(
    cbind(
      "2.5%" = intervalo.credibilidade[, 1],
      "Média" = sumario$Mean,
      "97.5%" = intervalo.credibilidade[, 2]
    ), row.names = row.names(intervalo.credibilidade)
  )
  
  return(resultados = list(sumario, tab))
}

      
# tabela_coeficientes -----------------------------------------------------------
      
tabela_coeficientes <- function(coef){
  plot_ly(
    type = 'table',
    columnorder = 1:ncol(coef),
    columnwidth = rep(80, ncol(coef)),
    header = list(
      values = coef %>% names %>% as.list,
      line = list(color = '#506784'),
      fill = list(color = "#1F8FFFB4"),
      align = rep('center',ncol(coef)),
      font = list(color = 'white', size = 15),
      height = 40
    ),
    cells = list(
      values = t(coef),
      line = list(color = '#506784'),
      fill = list(color = c("#1F8FFF58", rep('white',ncol(coef)-2),"#1F8FFF58")),
      align = rep('center',ncol(coef)),
      font = list(color = c('white',rep('#506784',ncol(coef)-2),'white'), size = 12),
      height = 30
    )) 
}

      
      
