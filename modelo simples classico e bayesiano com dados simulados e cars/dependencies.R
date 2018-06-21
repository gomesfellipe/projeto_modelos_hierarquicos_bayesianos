
# Cadeia ------------------------------------------------------------------

cadeia = function(df,name,p){
  par(mar = c(5,5,2,2),
      mfrow = c(ncol(df),1))
  for(i in 1:ncol(df)){
    plot(df[,i],type="l",ylab=name[i],xlab="Iterações", cex.lab=1.8,cex.main=2)
    abline(h=mean(p[i]), col="red",lwd=2)
    abline(h=quantile(df[,i],0.025),lty=3,col="blue",lwd=2)
    abline(h=quantile(df[,i],0.975),lty=3,col="blue",lwd=2)
  }
  par(mfrow=c(1,1))
}


# Histograma e densidade --------------------------------------------------

hist_den <- function(df, name = name,p){
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
  

# Autocorrelacao ----------------------------------------------------------

FAC <- function(df) {
  par(mar = c(5,5.5,2,2),
      mfrow = c(ncol(df),1))
  for(i in 1:ncol(df)){
    acf(x = df[,i],xlab="Defasagem",ylab="FAC",main="",
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
