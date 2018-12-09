source("serie_provision.r")



#########################################################################
#                              COMAUTO
#########################################################################

OthLiabData <- read.csv("comauto_pos.csv",header=TRUE, sep=",")

SumData <- ddply(OthLiabData,.(AccidentYear,DevelopmentYear,DevelopmentLag),summarise,IncurLoss=sum(IncurLoss_C-BulkLoss_C),CumPaidLoss=sum(CumPaidLoss_C), EarnedPremDIR=sum(EarnedPremDIR_C))
OL <- SumData[SumData$DevelopmentYear<1998,]
LossTri <- as.triangle(OL, origin="AccidentYear", dev="DevelopmentLag", value="CumPaidLoss")

prov_data<-function(data,n){
  Mdata<-matrix(data,nrow=n,ncol=n)
  prov=rep(0,n+1)
  temp=0
  for(i in 1:n){
    prov[i]=Mdata[i,n]-Mdata[i,n-i+1]
    
  }
  prov[n+1]=sum(prov[1:n])
  prov=matrix(prov,ncol=n+1)
  colnames(prov)<-c(1:n,"Total")
  
  return(prov)
}

#chain_ladder(LossTri)
comparaison(LossTri)
prov_data(SumData$CumPaidLoss,10)
points(1988:1997,prov_data(SumData$CumPaidLoss,10)[-11],type="b",lty=5,lwd=3,col="grey")

#-Provision globale compatible à ChainLadder
#-Provision annuelle plus proche du modèle multiplicatif

##########################################################################################
#                         Medmal
#########################################################################################
OthLiabData <- read.csv("medmal_pos.csv",header=TRUE, sep=",")

SumData <- ddply(OthLiabData,.(AccidentYear,DevelopmentYear,DevelopmentLag),summarise,IncurLoss=sum(IncurLoss_F2-BulkLoss_F2),CumPaidLoss=sum(CumPaidLoss_F2), EarnedPremDIR=sum(EarnedPremDIR_F2))
OL <- SumData[SumData$DevelopmentYear<1998,]
LossTri <- as.triangle(OL, origin="AccidentYear", dev="DevelopmentLag", value="CumPaidLoss")


comparaison(LossTri)
prov_data(SumData$CumPaidLoss,10)
points(1988:1997,prov_data(SumData$CumPaidLoss,10)[-11],type="b",lty=5,lwd=3,col="grey")

## ChainLadder plus propre

######################################################################################
#                GenInsLong avec chainLadder Fausse
######################################################################################


LossTri <- as.triangle(GenInsLong, origin="accyear", dev="devyear", value="incurred claims")
comparaison(LossTri)
prov_data(SumData$CumPaidLoss,10)
points(1988:1997,prov_data(SumData$CumPaidLoss,10)[-11],type="b",lty=5,lwd=3,col="grey")
#-- ChainLadder totallement abérrant ainsi que SARIMA.

######################################################################################
#     OTHER DATA
#####################################################################################

#comparaison(ABC)

#####################


comparaison(GenIns) 

