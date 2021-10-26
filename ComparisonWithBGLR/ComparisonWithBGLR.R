#Comparison between VIGoR and BGLR##############################################
library(VIGoR)
library(BGLR)
library(ROCR)

#Number of genotypes
N <- 1000

#Number of explanatory variables
##First variable
P1 <- c(1000,5000,10000,20000,50000)
##Second variable
P2 <- c(1000,5000,10000,20000,50000)

#Noise variance
Noise <- 0.25#Signal=1

#Number of repeats
Nsim <- 20


#Experiment 1###################################################################
#Proportions of explanatory variables with non-zero effects
Nzero1 <- 0.01
Nzero2 <- 0.001

#Objects for VIGoR
##Number of iterations
Ite.VIGoR.exp1 <- matrix(0, nr=Nsim, nc=length(P1))
##Calculation time
Sec.VIGoR.exp1 <- matrix(0, nr=Nsim, nc=length(P1))
##ACU
AUC.VIGoR.exp1 <- matrix(0, nr=Nsim, nc=length(P1))

#Objects for BGLR
##Calculation time
Sec.BGLR.exp1 <- matrix(0, nr=Nsim, nc=length(P1))
##AUC
AUC.BGLR.exp1 <- matrix(0, nr=Nsim, nc=length(P1))

#Run
for(i in 1:length(P1)){
  
  p1 <- P1[i]
  p2 <- P2[i]
  
  for(sim in 1:20){
    
    cat(i, sim, "\n")
    
    X1 <- matrix(sample(c(-1, 0, 1), N * p1, replace = TRUE), nr = N)
    X2 <- matrix(rnorm(N * p2), nr = N)
    
    B1 <- matrix(c(rep(1, p1*Nzero1), rep(0, p1-p1*Nzero1)), nc=1)
    B2 <- matrix(c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)), nc=1)
    
    Y <- as.vector(X1 %*% B1 + X2 %*% B2)
    E <- rnorm(N, 0, sqrt(var(Y) * Noise))
    Y <- Y + E
    
    X1.std <- scale(X1)
    X2.std <- scale(X2)
    
    #Calculate hyperparameter using hyperpara of VIGoR
    H1 <- hyperpara(X1.std, 0.5, "BayesC", Nzero1, "Var")
    H2 <- hyperpara(X2.std, 0.5, "BayesC", Nzero2, "Var")
    
    #RUn BGLR
    ETA <- list(list(model = "BayesC", 
                     X = X1.std, 
                     df0 = H1[1], S0 = H1[2], probIn = H1[3]),
                list(model = "BayesC", 
                     X = X2.std,
                     df0 = H2[1], S0 = H2[2], probIn = H2[3] ))
    
    start <- proc.time()[3]
    Result <- BGLR(y = Y, ETA = ETA, verbose = FALSE)
    end <- proc.time()[3]
    
    Sec.BGLR.exp1[sim, i] <- end - start
    
    ROCdata <- rbind(cbind(Result$ETA[[1]]$d, c(rep(1, p1*Nzero1), rep(0, p1-p1*Nzero1))),
                     cbind(Result$ETA[[2]]$d, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2))))
    ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
    Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
    AUC.temp <- performance(Pred, "auc")
    AUC.BGLR.exp1[sim, i] <- as.numeric(AUC.temp@y.values)    
    
    #Run VIGoR
    ETA <- list(list(model = "BayesC", 
                     X = X1.std, 
                     H = H1),
                list(model = "BayesC", 
                     X = X2.std,
                     H = H2))
    start <- proc.time()[3]
    Result <- vigor(Y, ETA, Verbose = FALSE, RandomIni = TRUE)
    end <- proc.time()[3]
    
    Ite.VIGoR.exp1[sim, i] <- length(Result$ResidualVar)
    Sec.VIGoR.exp1[sim, i] <- end - start
    
    ROCdata <- rbind(cbind(Result$ETA[[1]]$Rho, c(rep(1, p1*Nzero1), rep(0, p1-p1*Nzero1))),
                     cbind(Result$ETA[[2]]$Rho, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2))))
    ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
    Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
    AUC.temp <- performance(Pred, "auc")
    AUC.VIGoR.exp1[sim, i] <- as.numeric(AUC.temp@y.values)
  }
}


#Experiment 2###################################################################
Nzero1 <- 0.2
Nzero2 <- 0.001

Ite.VIGoR.exp2 <- matrix(0, nr=Nsim, nc=length(P1))
Sec.VIGoR.exp2 <- matrix(0, nr=Nsim, nc=length(P1))
AUC.VIGoR.exp2 <- matrix(0, nr=Nsim, nc=length(P1))

Sec.BGLR.exp2 <- matrix(0, nr=Nsim, nc=length(P1))
AUC.BGLR.exp2 <- matrix(0, nr=Nsim, nc=length(P1))

for(i in 1:length(P1)){
  
  p1 <- P1[i]
  p2 <- P2[i]
  
  for(sim in 1:20){
    
    cat(i, sim, "\n")
    
    X1 <- matrix(sample(c(-1, 0, 1), N * p1, replace = TRUE), nr = N)
    X2 <- matrix(rnorm(N * p2), nr = N)
    
    B1 <- matrix(c(rep(1, p1*Nzero1), rep(0, p1-p1*Nzero1)), nc=1)
    B2 <- matrix(c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)), nc=1)
    
    Y <- as.vector(X1 %*% B1 + X2 %*% B2)
    E <- rnorm(N, 0, sqrt(var(Y) * Noise))
    Y <- Y + E
    
    X1.std <- scale(X1)
    X2.std <- scale(X2)
    
    H1 <- hyperpara(X1.std, 0.5, "BL", Nzero1, "Var")
    H2 <- hyperpara(X2.std, 0.5, "BayesC", Nzero2, "Var")
    
    #Run BGLR
    ETA <- list(list(model = "BL", 
                     X = X1.std, 
                     shape = H1[1], rate = H1[2]),
                list(model = "BayesC", 
                     X = X2.std,
                     df0 = H2[1], S0 = H2[2], probIn = H2[3] ))
    
    start <- proc.time()[3]
    Result <- BGLR(y = Y, ETA = ETA, verbose = FALSE)
    end <- proc.time()[3]
    
    Sec.BGLR.exp2[sim, i] <- end - start
    
    ROCdata <- rbind(cbind(abs(Result$ETA[[1]]$b), c(rep(1, p1*Nzero1), rep(0, p1-p1*Nzero1))),
                     cbind(Result$ETA[[2]]$d, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2))))
    ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
    Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
    AUC.temp <- performance(Pred, "auc")
    AUC.BGLR.exp2[sim, i] <- as.numeric(AUC.temp@y.values)    
    
    #RUn VIGoR
    ETA <- list(list(model = "BL", 
                     X = X1.std, 
                     H = H1),
                list(model = "BayesC", 
                     X = X2.std,
                     H = H2))
    start <- proc.time()[3]
    Result <- vigor(Y, ETA, Verbose = FALSE, RandomIni = TRUE)
    end <- proc.time()[3]
    
    Ite.VIGoR.exp2[sim, i] <- length(Result$ResidualVar)
    Sec.VIGoR.exp2[sim, i] <- end - start
    
    ROCdata <- rbind(cbind(abs(Result$ETA[[1]]$Beta), c(rep(1, p1*Nzero1), rep(0, p1-p1*Nzero1))),
                     cbind(Result$ETA[[2]]$Rho, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2))))
    ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
    Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
    AUC.temp <- performance(Pred, "auc")
    AUC.VIGoR.exp2[sim, i] <- as.numeric(AUC.temp@y.values)
  }
}
rm(X1, X2, B1, B2, Y, E, X1.std, X2.std, ROCdata, Pred, AUC.temp)


#Plot###########################################################################
tiff("Result.tiff",width=17.5,height=17.5,unit="cm",res=600)
par(mfrow=c(2,2))
par(mar=c(3,3,2,0.5))
par(mgp=c(1.5,0.5,0))

#Experiment 1
s<-apply(Sec.BGLR.exp1,2,sd)
m<-apply(Sec.BGLR.exp1,2,mean)
plot(P1*2,m,type="o",xlab="",ylab="Sec.",main="Ex. 1",
     ylim=c(0,400))
arrows(P1*2,m,P1*2,m+s,length=0.1,angle=90)
arrows(P1*2,m,P1*2,m-s,length=0.1,angle=90)
s<-apply(Sec.VIGoR.exp1,2,sd)
m<-apply(Sec.VIGoR.exp1,2,mean)
points(P1*2,m,type="o",col=2,pch=8)
arrows(P1*2,m,P1*2,m+s,length=0.1,angle=90,col=2)
arrows(P1*2,m,P1*2,m-s,length=0.1,angle=90,col=2)

legend(0,400,col=c(2,1),pch=c(8,1),legend=c("VIGoR","BGLR"))

s<-apply(AUC.BGLR.exp1,2,sd)
m<-apply(AUC.BGLR.exp1,2,mean)
plot(P1*2,m,type="o",xlab="",ylab="AUC",main="Ex. 1",
     ylim=c(0.3,1))
arrows(P1*2,m,P1*2,m+s,length=0.1,angle=90)
arrows(P1*2,m,P1*2,m-s,length=0.1,angle=90)
s<-apply(AUC.VIGoR.exp1,2,sd)
m<-apply(AUC.VIGoR.exp1,2,mean)
points(P1*2,m,type="o",col=2,pch=8)
arrows(P1*2,m,P1*2,m+s,length=0.1,angle=90,col=2)
arrows(P1*2,m,P1*2,m-s,length=0.1,angle=90,col=2)

#Experiment 2
s<-apply(Sec.BGLR.exp2,2,sd)
m<-apply(Sec.BGLR.exp2,2,mean)
plot(P1*2,m,type="o",xlab="Number of explanatory variables",ylab="Sec.",main="Ex. 2",
     ylim=c(0,500))
arrows(P1*2,m,P1*2,m+s,length=0.1,angle=90)
arrows(P1*2,m,P1*2,m-s,length=0.1,angle=90)
s<-apply(Sec.VIGoR.exp2,2,sd)
m<-apply(Sec.VIGoR.exp2,2,mean)
points(P1*2,m,type="o",col=2,pch=8)
arrows(P1*2,m,P1*2,m+s,length=0.1,angle=90,col=2)
arrows(P1*2,m,P1*2,m-s,length=0.1,angle=90,col=2)

s<-apply(AUC.BGLR.exp2,2,sd)
m<-apply(AUC.BGLR.exp2,2,mean)
plot(P1*2,m,type="o",xlab="Number of explanatory variables",ylab="AUC",main="Ex. 2",
     ylim=c(0.4,1))
arrows(P1*2,m,P1*2,m+s,length=0.1,angle=90)
arrows(P1*2,m,P1*2,m-s,length=0.1,angle=90)
s<-apply(AUC.VIGoR.exp2,2,sd)
m<-apply(AUC.VIGoR.exp2,2,mean)
points(P1*2,m,type="o",col=2,pch=8)
arrows(P1*2,m,P1*2,m+s,length=0.1,angle=90,col=2)
arrows(P1*2,m,P1*2,m-s,length=0.1,angle=90,col=2)

dev.off()
