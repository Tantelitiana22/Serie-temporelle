require(TSA)
require(forecast)
library(ChainLadder)
library(plyr)
###################################################
#Liste de triangle de paiement cumulÃ©
#################################################
data(package="ChainLadder")
###################################################
# Fonction qui calcul la provision par utilisation
# de la projection diffÃ©renciÃ© dÃ©saisonnalisÃ©
####################################################
Provision_serie_df<-function(M,r=nrow(M),d_date){
  nl=nrow(M)
  Y_serie<-function(X,r=10){
    n=r*(1+r)/2
    M=X
    Y=rep(0,n)
    for(j in 1:r){
      I=((j-1)*j/2+1):((j-1)*(j+2)/2+1)
      tj=(j-1)*j/2+1
      for(t in I){
        Y[t]=M[abs(tj+j-t),abs(t-tj+1)]
      }
    }
    return(Y)
  }
  Y=Y_serie(M,r)
  for(i in 1:length(Y)){
    if(Y[i]<0){
    
      stop("Cette mÃ©thode ne fonctionne qu'avec des incrÃ©ments qui sont positifs!")
  }
  }
  Y=log(Y)
  require(TSA)
  require(forecast)
 
  dev.set(dev.prev())
  div=which.max(periodogram(Y)$spec)
  dev.off()
  floorD=floor(length(Y)/div)
  if(div==1){
    acf(Y,lag=50)
    stop("La sÃ©rie ne possÃ¨de pas de saisonnalitÃ©! Check ACF curve")
  }
  
  Z_t<-function(fit,L=90,r=10){
    Z=c()
    temp=0
    for(i in 1:L){
      if(i<=r+1){
        Z=c(Z,0)
      }else{
        temp=fit[i]-fit[i-r]-fit[i-1]+fit[i-r-1]
        Z=c(Z,temp)
      }
    }
    return(Z)
  }
  Z_reel<-Z_t(Y,L=nl*(nl+1)/2,r=floorD)
  
  
  p=suppressWarnings(pp.test(jitter(Z_reel)))
  if(p$p.value>0.05){
    return(list(Error="Les bruits ne sont pas stationnaire! Methode failed"))
  }
  arma_r<-auto.arima(Z_reel)
  ap=suppressWarnings(Box.test(arma_r$residuals))

  if(ap$p.value < 0.05){
    return(Error="Les bruits sont corrÃ©lÃ©s! Methode failed")
  }
  
  Prediction<-function(X,T,t,r=10,arma_mod){
    if(T>=t){
      return(X[t])
    }
    X_r=0
    temp=forecast(arma_mod,h=t-T)
    Z=sum(c(temp$mean,arma_mod$x[(r+2):T]))
    cste=X[11]-X[1]
    if(t-r>T){
      return(Z+cste+Prediction(X,T,t-r,r,arma_mod))
    }else{
      return(Z+cste+X[t-r])
    }
  }
  Val_pred<-function(X,T1,h,r_k=10,arma_m){
    res=c()
    for(i in 1:h){
      res=c(res,Prediction(X,T=T1,t=i+T,r=r_k,arma_mod=arma_m))
    }
    return(res)
  }
  pred=Val_pred(X=Y,T1=length(Y),r_k=floorD,h=nl*(nl-1),arma_m = arma_r)
  Y_k=c(Y,pred)
  
  M_serie<-function(Y,n,date=d_date){
    M=matrix(nrow=n,ncol=n)
    r=n
    for(j in 1:r){
      I=((j-1)*j/2+1):((j-1)*(j+2)/2+1)
      tj=(j-1)*j/2+1
      for(t in I){
        M[abs(tj+j-t),abs(t-tj+1)]=Y[t]
      }
    }
    
    for(k in 1:(n-1)){
      for(t in k:(n-1)){
        j=t+1
        i=n+k-t
        M[i,j]=Y[(n+1)*n/2+j+(k-1)*n]
      }
    } 
    rownames(M)<-1:n+date
    colnames(M)<-0:(n-1)
    return(M)
  }
  
  M_final=M_serie(exp(Y_k),nl)
  
  cum_res=M_final
  for(i in 1:(ncol(cum_res)-1)){
    cum_res[,i+1]=cum_res[,i+1]+cum_res[,i]
  }
  Prov_serie<-function(M,date=1988){
    nc<-ncol(M)
    nl<-nrow(M)
    prov=c()
    for(i in 1:nl){
      prov=c(prov,M[i,nc]-M[i,nl-i+1]) 
    }
    prov=matrix(prov,ncol=nl)
    prov=cbind(prov,sum(prov))
    colnames(prov)<-c(0:(nl-1)+date,"Total")
    return(prov)
  }
  prov=Prov_serie(cum_res,date=d_date+1)
  
  return(list(res_non_cum=M_final,res_cum=cum_res,provision=prov))
}

####################################################################
# Fonction permettant de remplir le tableau de paiement
# par les modÃ¨les de sÃ©rie additif. Si quadratique=TRUE
# On choisit alors de modÃ©liser avec un M.A Ã  tendance
# quadratique, sinon avec un M.A Ã  tendance linÃ©aire
####################################################################
Provision_serie_MA<-function(M,r=nrow(M),quadratique=FALSE,d_date=1987){
  nl=nrow(M)
  Y_serie<-function(X,r=10){
    n=r*(1+r)/2
    M=X
    Y=rep(0,n)
    for(j in 1:r){
      I=((j-1)*j/2+1):((j-1)*(j+2)/2+1)
      tj=(j-1)*j/2+1
      for(t in I){
        Y[t]=M[abs(tj+j-t),abs(t-tj+1)]
      }
    }
    return(Y)
  }
  Y=Y_serie(M,r)
  for(i in 1:length(Y)){
    if(Y[i]<0){
      stop("Cette fonction ne marche qu'avec des incrÃ©ments positifs")
    }
  }
  Y=log(Y)
  require(TSA)
  require(forecast)
  dev.set(dev.prev())
  div=which.max(periodogram(Y)$spec)
  dev.off()
  if(div==1){
    acf(Y,lag=50)
    stop("Reglement calendaire non saisonnaire! Check acf curve")
  }
  floorD=floor(length(Y)/div)
  
  
  res.lm=c()
  t1=1:length(Y)
  t2=t1^2
  a=0
  b=0
  c=0
  
  signe_difference<-function(vector){
    n=length(vector)
    U=rep(0,n-1) 
    if(n>1){
      if(vector[2]>vector[1]){
        U[1]=1 }
      for(i in 2:(n-1)){
        if(vector[i+1]>vector[i]){ 
          U[i]=1 
        }else{
          U[i]=0 }
      }
    }
    W=(sum(U)-(n-1)/4)/((n+1)/12)
    return(abs(W))
    
  }
  Moyenne_mobile<-function(X,t,Ta=12,da=4){
    if((t<=da%/%2)|(t>Ta-da%/%2)){
      return(0)
    }
    coeff=c()
    q=da%/%2
    
    #Cas oÃ¹ d est impaire
    if(da%%2==1){
      coeff=rep(1/da,2*q+1)
    }else{ # Sinon d est paire
      for(j in -q:q){
        if(abs(j)==q){
          coeff=c(coeff,1/(4*q))
        }else{
          coeff=c(coeff,1/(2*q))
        }
      }
    }
    res=rep(0,2*q+1)
    for(i in -q:q){
      res[i+q+1]=coeff[i+q+1]*X[t-i]
    }
    return(sum(res))
  }
  
  saison2_fontion<-function(vect,T=12,d=4){
    J=floor(T/d)-1
    X<-vect
    W=rep(0,d)
    S1=rep(0,d)
    # Application de la formule 
    for(k in 1:d){
      temp=0
      for(j in 0:J){
        temp=temp+X[k+d*j]-Moyenne_mobile(X,t=k+d*j,Ta=T,da=d)
      }
      W[k]=temp/J
    }
    
    S1=W-sum(W)/d
    # Condition qui nous permet d'avoir la somme des saisons S_k Ã©gale Ã  0.
    S1=S1-sum(S1)/d
    res=c()
    #On repete (T/d) fois la tendance.
    for(i in 1:(floor(T/d))){
      res=c(res,S1)
    }
    return(res)
  }
  Equation=rep(0,length(t1))
  if(signe_difference(Y)>1.96){
    if(quadratique==TRUE){
      res.lm=lm(Y~t1+t2)
      a=res.lm$coefficients[1]
      b=res.lm$coefficients[2]
      c=res.lm$coefficients[3]
      Equation=a+b*t1+c*t2
    }else{
      res.lm=lm(Y~t1)
      a=res.lm$coefficients[1]
      b=res.lm$coefficients[2]
      Equation=a+b*t1
    
    }
  }else{
    print(" ")
    warning(cat("Pas de tendance détécté:SD W=",signe_difference(Y),"<1.96"))
    print(" ")
  }
  Extraction=Y-Equation  
  saison=saison2_fontion(Extraction,T=length(Y),d=floorD)
  l_sais=length(saison)
  l_equa=length(Equation)
  len=0
  bruit=0
  
  if(l_equa>l_sais){
    len=l_equa-l_sais     
    bruit=Y-(Equation+c(saison,saison[1:len]))
  }else{
    bruit=Y-(Equation+saison)
  }

  p=suppressWarnings(pp.test(bruit))
  if(p$p.value>0.05){
    stop("Les bruits ne sont pas stationnaire! Methode failed")
  }
  arma<-auto.arima(bruit)
  at=suppressWarnings(Box.test(arma$residuals))
  if(at$p.value<0.05){
    stop("Les bruits sont corrélés! Methode failed")
  }

  h_pred=nl*(nl-1)  
  b_pred=forecast(arma,h=h_pred)
  h1=(nl*(nl+1)/2+1):(nl*(nl+1)/2+h_pred)
  h2=h1^2 
  Eq=a+b*h1+c*h2
  dsaison=c()
  dsaison=c(dsaison,saison[(len+1):l_sais])
  while(length(dsaison)<h_pred){
    if(h_pred>length(dsaison)+l_sais){
      dsaison=c(dsaison,saison[1:l_sais])
    }else{
      dsaison=c(dsaison,saison[1:(h_pred-length(dsaison))])
    }
  }

  Y_predict=Eq+dsaison+b_pred$mean
  Y_fin=c(Y,Y_predict)

  M_serie<-function(Y,n,date=d_date){
    M=matrix(nrow=n,ncol=n)
    r=n
    for(j in 1:r){
      I=((j-1)*j/2+1):((j-1)*(j+2)/2+1)
      tj=(j-1)*j/2+1
      for(t in I){
        M[abs(tj+j-t),abs(t-tj+1)]=Y[t]
      }
    }
    
    for(k in 1:(n-1)){
      for(t in k:(n-1)){
        j=t+1
        i=n+k-t
        M[i,j]=Y[(n+1)*n/2+j+(k-1)*n]
      }
    } 
    rownames(M)<-1:n+date
    colnames(M)<-0:(n-1)
    return(M)
  } 
  M_final=M_serie(exp(Y_fin),nl)
  
  cum_res=M_final
  for(i in 1:(ncol(cum_res)-1)){
    cum_res[,i+1]=cum_res[,i+1]+cum_res[,i]
  }
  Prov_serie<-function(M,date=1988){
    nc<-ncol(M)
    nl<-nrow(M)
    prov=c()
    for(i in 1:nl){
      prov=c(prov,M[i,nc]-M[i,nl-i+1]) 
    }
    prov=matrix(prov,ncol=nl)
    prov=cbind(prov,sum(prov))
    colnames(prov)<-c(0:(nl-1)+date,"Total")
    return(prov)
  }
  prov=Prov_serie(cum_res,date=d_date+1)
  
  return(list(res_non_cum=M_final,res_cum=cum_res,provision=prov))
  
}

##################################################################
# MÃ©thode de Chain Ladder
#################################################################
chain_ladder<-function(M){
  f=rep(0,ncol(M)-1)
  M[is.na(M)] <- 0
  for(i in 1:(ncol(M)-1)){
    f[i]=sum(M[1:(nrow(M)-i),i+1])/sum(M[1:(nrow(M)-i),i]) 
  }
  triangle=M
  for( i in 1:(ncol(M)-1)){
    triangle[(nrow(triangle)-i+1):(nrow(triangle)),i+1]=
      f[i]*triangle[(nrow(triangle)-i+1):(nrow(triangle)),i]
  }
  Non_cum=triangle
  Non_cum[,2:ncol(Non_cum)]=Non_cum[,2:ncol(Non_cum)]-Non_cum[,1:(ncol(Non_cum)-1)]
  
  prov=rep(0,nrow(Non_cum)+1)
  for(j in 1:nrow(Non_cum)){
    if(j==1){
      prov[j]=0
    }else{
      prov[j]=sum(Non_cum[j,(ncol(Non_cum)-j+2):ncol(Non_cum)])
    }
  }
  prov[nrow(Non_cum)+1]=sum(prov)
  prov=matrix(prov,ncol=ncol(Non_cum)+1,nrow=1)
  colnames(prov)<-c(rownames(Non_cum),"Total")
  return(list(facteur=f,PAID=floor(triangle),non_cumul=floor(Non_cum)
              ,provision=floor(prov)))
}

######################################################################
#MÃ©thode de London Chain Ladder
#####################################################################

London_chain_ladder<-function(M){
  f=rep(0,ncol(M)-1)
  a=rep(0,ncol(M)-1)
  M[is.na(M)] <- 0
  temp=0
  for(i in 1:(ncol(M)-1)){
    if(i+1!=ncol(M)){
      temp=lm(M[1:(nrow(M)-i),i+1]~M[1:(nrow(M)-i),i])
      a[i]=coef(temp)[1]
      f[i]=coef(temp)[2]
    }else{
      a[i]=0
      f[i]=sum(M[1:(nrow(M)-i),i+1])/sum(M[1:(nrow(M)-i),i])
    }
  }
  triangle=M
  for( i in 1:(ncol(M)-1)){
    triangle[(nrow(triangle)-i+1):(nrow(triangle)),i+1]=
      f[i]*triangle[(nrow(triangle)-i+1):(nrow(triangle)),i]+a[i]
  }
  Non_cum=triangle
  Non_cum[,2:ncol(Non_cum)]=Non_cum[,2:ncol(Non_cum)]-Non_cum[,1:(ncol(Non_cum)-1)]
  
  prov=rep(0,nrow(Non_cum)+1)
  for(j in 1:nrow(Non_cum)){
    if(j==1){
      prov[j]=0
    }else{
      prov[j]=sum(Non_cum[j,(ncol(Non_cum)-j+2):ncol(Non_cum)])
    }
  }
  prov[nrow(Non_cum)+1]=sum(prov)
  prov=matrix(prov,ncol=ncol(Non_cum)+1,nrow=1)
  colnames(prov)<-c(rownames(Non_cum),"Total")
  
  return(list(facteur=f,intercept=a,PAID=floor(triangle),non_cumul=floor(Non_cum)
              ,provision=floor(prov)))
}


comparaison<-function(M,d_date=1987){

 print("provison avec Chain Ladder")
 print(chain_ladder(M)$provision)
 print("provison avec London Chain Ladder")
 print(London_chain_ladder(M)$provision)
 non_cum<-M
 non_cum[,2:ncol(non_cum)]<- -non_cum[,1:(ncol(non_cum)-1)]+non_cum[,2:ncol(non_cum)]

 print("provison avec  SARIMA")
 print(floor(Provision_serie_df(non_cum,d_date = 1987)$provision))
 print("provison avec  M.A  tendance linéaire")
 print(floor(Provision_serie_MA(non_cum)$provision))
 print("provison avec  M.A tendance quadratique")
 print(floor(Provision_serie_MA(non_cum,quadratique = TRUE)$provision))

 plot(1:ncol(M)+d_date,London_chain_ladder(M)$provision[-(nrow(M)+1)],type="b",ylab="Provision",xlab="Date",
      lwd=4,col=4,main="Provisionnement avec les différentes méthodes")
 points(1:ncol(M)+d_date,chain_ladder(M)$provision[-(nrow(M)+1)],type="b",lwd="4",col=5)
 points(1:ncol(M)+d_date,Provision_serie_df(non_cum,d_date = 1987)$provision[-(nrow(M)+1)],type="b",lwd=7,lty=8,col=6)
 points(1:ncol(M)+d_date,Provision_serie_MA(non_cum)$provision[-(nrow(M)+1)],type="b",lwd="4",col=2)
 points(1:ncol(M)+d_date,Provision_serie_MA(non_cum,quadratique = TRUE)$provision[-(nrow(M)+1)],type="b",lwd=3,col=9)
 
 grid (10,10, lty = 6, col = "cornsilk2")
 legend("topleft",c("Chain Ladder","London Chain Ladder",
                    "Projection SARIMA","M.A linéaire","M.A quadratique"),
        col=c(4,5,6,2,9),lty=c(1,1,1,1),lwd=c(2.5,2.5,2.5,2.5),cex=0.8)
}





signe_difference<-function(vector){
  n=length(vector)
  U=rep(0,n-1) 
  if(n>1){
    if(vector[2]>vector[1]){
      U[1]=1 }
  for(i in 2:(n-1)){
    if(vector[i+1]>vector[i]){ 
      U[i]=1 
    }else{
      U[i]=0 }
  }
  }
  W=(sum(U)-(n-1)/4)/((n+1)/12)
  res="" 
  if(abs(W)>1.96){
    return(list(W_stat=abs(W),commentaire="La série possède une tendance")) 
  }else{ 
      return(list(W_stat=abs(W),commentaire="La série ne possède pas de tendance"))
  }

}

