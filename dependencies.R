
# Pacotes -----------------------------------------------------------------

packages = c('dplyr', 'ggplot2','gridExtra', 'dplyr','purrr', 'xtable','stringr', 'ggExtra','plotly')

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    suppressMessages(library(package, character.only=T))
  }
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
      font = list(color = c('#506784'), size = 12),
      height = 30
    )) 
}

      
      
