#Packages#######################################################################
library(VIGoR)
library(BGLR)
library(ROCR)
library(snow)
#This script conducts parallel computation using snow, which requires large memory
#in particular when the number of explanatory variables is large


#Simulation conditions##########################################################
#Number of samples
N <- 1000

#Number of explanatory variables
P1 <- c(1000,5000,10000,20000,50000)
P2 <- c(1000,5000,10000,20000,50000)

#Residual variance
Noise <- 0.25

#Number of replications
Nsim <- 20


#Exp. 1########################################################################
#Proportion of variables with non-zero effects
Nzero1 <- 0.01
Nzero2 <- 0.001

#Create directories
for(sim in 1:Nsim) dir.create(as.character(sim))

#Create data files in each directory
for(i in 1:1){
  
  p1 <- P1[i]
  p2 <- P2[i]
  
  for(sim in 1:20){
    
    X1 <- scale(matrix(rnorm(2 * N * p1), nr = 2 * N))
    X2 <- scale(matrix(rnorm(2 * N * p2), nr = 2 * N))
    
    B1 <- matrix(c(rnorm(p1*Nzero1), rep(0, p1-p1*Nzero1)), nc=1)
    B2 <- matrix(c(rnorm(p2*Nzero2), rep(0, p2-p2*Nzero2)), nc=1)
    
    U1 <- scale(as.numeric(X1 %*% B1))
    U2 <- scale(as.numeric(X2 %*% B2))
    U <- U1 + U2
    E <- rnorm(N, 0, sqrt(var(U) * Noise))
    Y <- U + E
    
    write.csv(cbind(X1, X2), paste(sim, "/Exp1_", i, "_", sim, "_X.csv", sep = ""), row.names = FALSE)
    write.csv(cbind(B1, B2), paste(sim, "/Exp1_", i, "_", sim, "_B.csv", sep = ""), row.names = FALSE)
    write.csv(cbind(Y, U), paste(sim, "/Exp1_", i, "_", sim, "_Y.csv", sep = ""), row.names = FALSE)
  }
}

#Run
##Parallel computation using snow
##20 replications are run using 5 clusters.
Stats.Exp1.BayesCBayesC <- NULL
cl <- makeCluster(5, type = "SOCK")
clusterExport(cl, c("BGLR", "vigor", "hyperpara", "predict_vigor", "prediction", "performance",
                    "N", "Nzero1", "Nzero2", "Noise"))

for(i in 1:length(P1)){
  
  p1 <- P1[i]
  p2 <- P2[i]
  clusterExport(cl, c("p1", "p2", "i"))
  
  for(block in 1:4){
    
    cat(i, p1, p2, block, "\n")
    
    Data.list <- as.list(numeric(5))
    for(j in 1:5){
      d <- (block - 1) * 5 + j
      X <- matrix(scan(paste(d, "/Exp1_", i, "_", d, "_X.csv", sep = ""), skip = 1, sep = ",", quiet = TRUE),
                  nrow = 2000, byrow = TRUE)
      Y <- as.matrix(read.csv(paste(d, "/Exp1_", i, "_", d, "_Y.csv", sep = "")))
      Data.list[[j]] <- list(d = d,
                             X1 = X[, 1:p1],
                             X2 = X[, (p1 + 1):(p1 + p2)],
                             U = Y[, 2],
                             Y = Y[, 1])
    }
    
    Result <- parLapply(cl,
                        Data.list,
                        function(Data){
                          
                          Result.sim1 <- numeric(39)
                          
                          #Hyperparameters
                          H1 <- hyperpara(Data$X1[1:N, ], 1.0, "BayesC", Nzero1, "Var")
                          H2 <- hyperpara(Data$X2[1:N, ], 1.0, "BayesC", Nzero2, "Var")
                          
                          ##BGLR
                          ETA <- list(list(model = "BayesC",
                                           X = Data$X1[1:N, ],
                                           df0 = H2[1], S0 = H2[2], probIn = H2[3], counts = 1e+6),
                                      list(model = "BayesC",
                                           X = Data$X2[1:N, ],
                                           df0 = H2[1], S0 = H2[2], probIn = H2[3], counts = 1e+6))
                          
                          ###1500
                          start <- proc.time()[3]
                          Result <- BGLR(y = Data$Y[1:N], ETA = ETA, verbose = FALSE, burnIn = 500, nIter = 1500, thin = 5,
                                         saveAt = paste(Data$d, "/", sep = ""))
                          end <- proc.time()[3]
                          
                          Result.sim1[1] <- end - start
                          
                          ####AUC
                          ROCdata <- cbind(Result$ETA[[1]]$d, c(rep(1, p1*Nzero1), rep(0, p1-p1*Nzero1)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim1[2] <- as.numeric(AUC.temp@y.values)
                          
                          ROCdata <- cbind(Result$ETA[[2]]$d, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim1[3] <- as.numeric(AUC.temp@y.values)
                          
                          ####Rhat
                          s <- scan(paste(Data$d, "varE.dat", sep = "/"))
                          s <- matrix(s[(length(s)-199):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim1[4] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- read.table(paste(Data$d, "ETA_1_parBayesC.dat", sep = "/"))[, 1]
                          s <- matrix(s[(length(s)-199):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim1[5] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- read.table(paste(Data$d, "ETA_2_parBayesC.dat", sep = "/"))[, 1]
                          s <- matrix(s[(length(s)-199):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim1[6] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- scan(paste(Data$d, "mu.dat", sep = "/"))
                          s <- matrix(s[(length(s)-199):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim1[7] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          ####Prediction
                          temp <- Result$mu + 
                            colSums(t(Data$X1[(N + 1):(2 * N), ]) * Result$ETA[[1]]$b) + 
                            colSums(t(Data$X2[(N + 1):(2 * N), ]) * Result$ETA[[2]]$b)
                          Result.sim1[8] <- cor(Data$U[(N + 1):(2 * N)], temp)
                          
                          ###15000
                          start <- proc.time()[3]
                          Result <- BGLR(y = Data$Y[1:N], ETA = ETA, verbose = FALSE, burnIn = 5000, nIter = 15000, thin = 10,
                                         saveAt = paste(Data$d, "/", sep = ""))
                          end <- proc.time()[3]
                          
                          Result.sim1[9] <- end - start
                          
                          ####AUC
                          ROCdata <- cbind(Result$ETA[[1]]$d, c(rep(1, p1*Nzero1), rep(0, p1-p1*Nzero1)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim1[10] <- as.numeric(AUC.temp@y.values)
                          
                          ROCdata <- cbind(Result$ETA[[2]]$d, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim1[11] <- as.numeric(AUC.temp@y.values)
                          
                          ####Rhat
                          s <- scan(paste(Data$d, "varE.dat", sep = "/"))
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim1[12] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- read.table(paste(Data$d, "ETA_1_parBayesC.dat", sep = "/"))[, 1]
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim1[13] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- read.table(paste(Data$d, "ETA_2_parBayesC.dat", sep = "/"))[, 1]
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim1[14] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- scan(paste(Data$d, "mu.dat", sep = "/"))
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim1[15] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          ####Prediction
                          temp <- Result$mu + 
                            colSums(t(Data$X1[(N + 1):(2 * N), ]) * Result$ETA[[1]]$b) + 
                            colSums(t(Data$X2[(N + 1):(2 * N), ]) * Result$ETA[[2]]$b)
                          Result.sim1[16] <- cor(Data$U[(N + 1):(2 * N)], temp)
                          
                          ###30000
                          start <- proc.time()[3]
                          Result <- BGLR(y = Data$Y[1:N], ETA = ETA, verbose = FALSE, burnIn = 20000, nIter = 30000, thin = 10,
                                         saveAt = paste(Data$d, "/", sep = ""))
                          end <- proc.time()[3]
                          
                          Result.sim1[17] <- end - start
                          
                          ####AUC
                          ROCdata <- cbind(Result$ETA[[1]]$d, c(rep(1, p1*Nzero1), rep(0, p1-p1*Nzero1)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim1[18] <- as.numeric(AUC.temp@y.values)
                          
                          ROCdata <- cbind(Result$ETA[[2]]$d, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim1[19] <- as.numeric(AUC.temp@y.values)
                          
                          ####Rhat
                          s <- scan(paste(Data$d, "varE.dat", sep = "/"))
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim1[20] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- read.table(paste(Data$d, "ETA_1_parBayesC.dat", sep = "/"))[, 1]
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim1[21] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- read.table(paste(Data$d, "ETA_2_parBayesC.dat", sep = "/"))[, 1]
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim1[22] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- scan(paste(Data$d, "mu.dat", sep = "/"))
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim1[23] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          ####Prediction
                          temp <- Result$mu + 
                            colSums(t(Data$X1[(N + 1):(2 * N), ]) * Result$ETA[[1]]$b) + 
                            colSums(t(Data$X2[(N + 1):(2 * N), ]) * Result$ETA[[2]]$b)
                          Result.sim1[24] <- cor(Data$U[(N + 1):(2 * N)], temp)
                          
                          ##VIGoR
                          ETA <- list(list(model = "BayesC",
                                           X = Data$X1[1:N, ],
                                           H = H1),
                                      list(model = "BayesC",
                                           X = Data$X2[1:N, ],
                                           H = H2))
                          
                          ###1e-4
                          start <- proc.time()[3]
                          Result <- vigor(Data$Y[1:N], ETA, Verbose = FALSE, RandomIni = TRUE, Thresholdvalue = 1e-4, Maxiteration = 10000)
                          end <- proc.time()[3]
                          
                          Result.sim1[25] <- end - start
                          Result.sim1[26] <- length(Result$ResidualVar)
                          
                          ####AUC
                          ROCdata <- cbind(Result$ETA[[1]]$Rho, c(rep(1, p1*Nzero1), rep(0, p1-p1*Nzero1)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim1[27] <- as.numeric(AUC.temp@y.values)
                          
                          ROCdata <- cbind(Result$ETA[[2]]$Rho, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim1[28] <- as.numeric(AUC.temp@y.values)
                          
                          ####Prediction
                          temp <- predict_vigor(Result, list(Data$X1[(N + 1):(2 * N), ], Data$X2[(N + 1):(2 * N), ]))
                          Result.sim1[29] <- cor(Data$U[(N + 1):(2 * N)], temp)
                          
                          ###1e-5
                          start <- proc.time()[3]
                          Result <- vigor(Data$Y[1:N], ETA, Verbose = FALSE, RandomIni = TRUE, Thresholdvalue = 1e-5, Maxiteration = 10000)
                          end <- proc.time()[3]
                          
                          Result.sim1[30] <- end - start
                          Result.sim1[31] <- length(Result$ResidualVar)
                          
                          ####AUC
                          ROCdata <- cbind(Result$ETA[[1]]$Rho, c(rep(1, p1*Nzero1), rep(0, p1-p1*Nzero1)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim1[32] <- as.numeric(AUC.temp@y.values)
                          
                          ROCdata <- cbind(Result$ETA[[2]]$Rho, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim1[33] <- as.numeric(AUC.temp@y.values)
                          
                          ####Prediction
                          temp <- predict_vigor(Result, list(Data$X1[(N + 1):(2 * N), ], Data$X2[(N + 1):(2 * N), ]))
                          Result.sim1[34] <- cor(Data$U[(N + 1):(2 * N)], temp)
                          
                          ###1e-6
                          start <- proc.time()[3]
                          Result <- vigor(Data$Y[1:N], ETA, Verbose = FALSE, RandomIni = TRUE, Thresholdvalue = 1e-6, Maxiteration = 10000)
                          end <- proc.time()[3]
                          
                          Result.sim1[35] <- end - start
                          Result.sim1[36] <- length(Result$ResidualVar)
                          
                          ####AUC
                          ROCdata <- cbind(Result$ETA[[1]]$Rho, c(rep(1, p1*Nzero1), rep(0, p1-p1*Nzero1)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim1[37] <- as.numeric(AUC.temp@y.values)
                          
                          ROCdata <- cbind(Result$ETA[[2]]$Rho, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim1[38] <- as.numeric(AUC.temp@y.values)
                          
                          ####Prediction
                          temp <- predict_vigor(Result, list(Data$X1[(N + 1):(2 * N), ], Data$X2[(N + 1):(2 * N), ]))
                          Result.sim1[39] <- cor(Data$U[(N + 1):(2 * N)], temp)
                          
                          Result.sim1
                        }#function
    )#parLapply
    
    temp <- matrix(unlist(Result), nr = 39)
    Stats.Exp1.BayesCBayesC <- rbind(Stats.Exp1.BayesCBayesC, temp)
    
  }#block
}#i
stopCluster(cl)

#See the averages
for(i in 1:length(P1)){
  
  temp <- Stats.Exp1.BayesCBayesC[((i - 1) * 39 * 4 + 1):(i * 39 * 4),]
  temp <- apply(matrix(as.numeric(temp), nr = 39), 1, mean)
  
  cat("#BGLR","\n")
  cat("#", P1[i] * 2, temp[1:8], "\n")#1500
  cat("#", P1[i] * 2, temp[9:16], "\n")#20000
  cat("#", P1[i] * 2, temp[17:24], "\n")#30000
  cat("#VIGoR","\n")  
  cat("#", P1[i] * 2, temp[25:29], "\n")#1e-4
  cat("#", P1[i] * 2, temp[30:34], "\n")#1e-5
  cat("#", P1[i] * 2, temp[35:39], "\n\n")#1e-6
}


#Exp 2###########################################################################
Nzero1 <- 1
Nzero2 <- 0.001

for(i in 1:length(P1)){
  
  p1 <- P1[i]
  p2 <- P2[i]
  
  for(sim in 1:20){
    
    X1 <- scale(matrix(sample(c(-1, 0, 1), 2 * N * p1, replace = TRUE, prob = c(0.25, 0.5, 0.25)), nr = 2 * N))
    X2 <- scale(matrix(rnorm(2 * N * p2), nr = 2 * N))
    
    B1 <- matrix(c(rnorm(p1 * Nzero1), rep(0, p1 - p1 * Nzero1)), nc = 1)
    B2 <- matrix(c(rnorm(p2 * Nzero2), rep(0, p2 - p2 * Nzero2)), nc = 1)
    
    U1 <- scale(as.numeric(X1 %*% B1))
    U2 <- scale(as.numeric(X2 %*% B2))
    U <- U1 + U2
    E <- rnorm(N, 0, sqrt(var(U) * Noise))
    Y <- U + E
    
    write.csv(cbind(X1, X2), paste(sim, "/Exp2_", i, "_", sim, "_X.csv", sep = ""), row.names = FALSE)
    write.csv(cbind(B1, B2), paste(sim, "/Exp2_", i, "_", sim, "_B.csv", sep = ""), row.names = FALSE)
    write.csv(cbind(Y, U), paste(sim, "/Exp2_", i, "_", sim, "_Y.csv", sep = ""), row.names = FALSE)
  }
}

Stats.Exp2.BRRBayesB <- NULL
cl <- makeCluster(5, type = "SOCK")
clusterExport(cl, c("BGLR", "vigor", "hyperpara", "predict_vigor", "prediction", "performance",
                    "N", "Nzero1", "Nzero2", "Noise"))

for(i in 1:length(P1)){
  
  p1 <- P1[i]
  p2 <- P2[i]
  clusterExport(cl, c("p1", "p2", "i"))
  
  for(block in 1:4){
    
    cat(i, p1, p2, block, "\n")
    
    Data.list <- as.list(numeric(5))
    for(j in 1:5){
      d <- (block - 1) * 5 + j
      X <- matrix(scan(paste(d, "/Exp2_", i, "_", d, "_X.csv", sep = ""), skip = 1, sep = ",", quiet = TRUE),
                  nrow = 2000, byrow = TRUE)
      Y <- as.matrix(read.csv(paste(d, "/Exp2_", i, "_", d, "_Y.csv", sep = "")))
      B <- as.matrix(read.csv(paste(d, "/Exp2_", i, "_", d, "_B.csv", sep = "")))
      Data.list[[j]] <- list(d = d,
                             X1 = X[, 1:p1],
                             X2 = X[, (p1 + 1):(p1 + p2)],
                             B1 = B[, 1],
                             U = Y[, 2],
                             Y = Y[, 1])
    }
    
    Result <- parLapply(cl,
                        Data.list,
                        function(Data){
                          
                          Result.sim2 <- numeric(39)
                          
                          H1 <- hyperpara(Data$X1[1:N, ], 1.0, "BRR", Nzero1, "Var")
                          H2 <- hyperpara(Data$X2[1:N, ], 1.0, "BayesB", Nzero2, "Var")
                          
                          ##BGLR
                          ETA <- list(list(model = "BRR",
                                           X = Data$X1[1:N, ],
                                           df0 = H1[1], S0 = H1[2]),
                                      list(model = "BayesB",
                                           X = Data$X2[1:N, ],
                                           df0 = H2[1], S0 = H2[2], probIn = H2[3], counts = 1e+6))
                          
                          ###1500
                          start <- proc.time()[3]
                          Result <- BGLR(y = Data$Y[1:N], ETA = ETA, verbose = FALSE, burnIn = 500, nIter = 1500, thin = 5,
                                         saveAt = paste(Data$d, "/", sep = ""))
                          end <- proc.time()[3]
                          
                          Result.sim2[1] <- end - start
                          
                          ####Estimation accuracy
                          Result.sim2[2] <- cor(Result$ETA[[1]]$b, as.numeric(Data$B1))
                          
                          ####AUC
                          ROCdata <- cbind(Result$ETA[[2]]$d, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim2[3] <- as.numeric(AUC.temp@y.values)
                          
                          ####Rhat
                          s <- scan(paste(Data$d, "varE.dat", sep = "/"))
                          s <- matrix(s[(length(s)-199):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim2[4] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- read.table(paste(Data$d, "ETA_1_varB.dat", sep = "/"))[, 1]
                          s <- matrix(s[(length(s)-199):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim2[5] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- scan(paste(Data$d, "mu.dat", sep = "/"))
                          s <- matrix(s[(length(s)-199):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim2[7] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          ####Prediction accuracy
                          temp <- Result$mu + 
                            colSums(t(Data$X1[(N + 1):(2 * N), ]) * Result$ETA[[1]]$b) + 
                            colSums(t(Data$X2[(N + 1):(2 * N), ]) * Result$ETA[[2]]$b)
                          Result.sim2[8] <- cor(Data$U[(N + 1):(2 * N)], temp)
                          
                          ###15000
                          start <- proc.time()[3]
                          Result <- BGLR(y = Data$Y[1:N], ETA = ETA, verbose = FALSE, burnIn = 5000, nIter = 15000, thin = 10,
                                         saveAt = paste(Data$d, "/", sep = ""))
                          end <- proc.time()[3]
                          
                          Result.sim2[9] <- end - start
                          
                          ####Estimation accuracy
                          Result.sim2[10] <- cor(Result$ETA[[1]]$b, as.numeric(Data$B1))
                          
                          ####AUC
                          ROCdata <- cbind(Result$ETA[[2]]$d, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim2[11] <- as.numeric(AUC.temp@y.values)
                          
                          ####Rhat
                          s <- scan(paste(Data$d, "varE.dat", sep = "/"))
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim2[12] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- read.table(paste(Data$d, "ETA_1_varB.dat", sep = "/"))[, 1]
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim2[13] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- scan(paste(Data$d, "mu.dat", sep = "/"))
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim2[15] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          ####Prediction accuracy
                          temp <- Result$mu + 
                            colSums(t(Data$X1[(N + 1):(2 * N), ]) * Result$ETA[[1]]$b) + 
                            colSums(t(Data$X2[(N + 1):(2 * N), ]) * Result$ETA[[2]]$b)
                          Result.sim2[16] <- cor(Data$U[(N + 1):(2 * N)], temp)
                          
                          ###30000
                          start <- proc.time()[3]
                          Result <- BGLR(y = Data$Y[1:N], ETA = ETA, verbose = FALSE, burnIn = 20000, nIter = 30000, thin = 10,
                                         saveAt = paste(Data$d, "/", sep = ""))
                          end <- proc.time()[3]
                          
                          Result.sim2[17] <- end - start
                          
                          ####Estimation accuracy
                          Result.sim2[18] <- cor(Result$ETA[[1]]$b, as.numeric(Data$B1))
                          
                          ####AUC
                          ROCdata <- cbind(Result$ETA[[2]]$d, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim2[19] <- as.numeric(AUC.temp@y.values)
                          
                          ####Rhat
                          s <- scan(paste(Data$d, "varE.dat", sep = "/"))
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim2[20] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- read.table(paste(Data$d, "ETA_1_varB.dat", sep = "/"))[, 1]
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim2[21] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          s <- scan(paste(Data$d, "mu.dat", sep = "/"))
                          s <- matrix(s[(length(s)-999):length(s)], nc=2)
                          B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                          W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                          Result.sim2[23] <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                          
                          ####Prediction accuracy
                          temp <- Result$mu + 
                            colSums(t(Data$X1[(N + 1):(2 * N), ]) * Result$ETA[[1]]$b) + 
                            colSums(t(Data$X2[(N + 1):(2 * N), ]) * Result$ETA[[2]]$b)
                          Result.sim2[24] <- cor(Data$U[(N + 1):(2 * N)], temp)
                          
                          ##VIGoR
                          ETA <- list(list(model = "BRR",
                                           X = Data$X1[1:N, ],
                                           H = H1),
                                      list(model = "BayesB",
                                           X = Data$X2[1:N, ],
                                           H = H2))

                          ###1e-4
                          start <- proc.time()[3]
                          Result <- vigor(Data$Y[1:N], ETA, Verbose = FALSE, RandomIni = TRUE, Thresholdvalue = 1e-4, Maxiteration = 10000)
                          end <- proc.time()[3]
                          
                          Result.sim2[25] <- end - start
                          Result.sim2[26] <- length(Result$ResidualVar)
                          
                          ####Estimation accuracy
                          Result.sim2[27] <- cor(Result$ETA[[1]]$Beta, as.numeric(Data$B1))
                          
                          ####AUC
                          ROCdata <- cbind(Result$ETA[[2]]$Rho, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim2[28] <- as.numeric(AUC.temp@y.values)
                          
                          ####Prediction accuracy
                          temp <- predict_vigor(Result, list(Data$X1[(N + 1):(2 * N), ], Data$X2[(N + 1):(2 * N), ]))
                          Result.sim2[29] <- cor(Data$U[(N + 1):(2 * N)], temp)
                          
                          ###1e-5
                          start <- proc.time()[3]
                          Result <- vigor(Data$Y[1:N], ETA, Verbose = FALSE, RandomIni = TRUE, Thresholdvalue = 1e-5, Maxiteration = 10000)
                          end <- proc.time()[3]
                          
                          Result.sim2[30] <- end - start
                          Result.sim2[31] <- length(Result$ResidualVar)
                          
                          ####Estimation accuracy
                          Result.sim2[32] <- cor(Result$ETA[[1]]$Beta, as.numeric(Data$B1))
                          
                          ####AUC
                          ROCdata <- cbind(Result$ETA[[2]]$Rho, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim2[33] <- as.numeric(AUC.temp@y.values)
                          
                          ####Prediction accuracy
                          temp <- predict_vigor(Result, list(Data$X1[(N + 1):(2 * N), ], Data$X2[(N + 1):(2 * N), ]))
                          Result.sim2[34] <- cor(Data$U[(N + 1):(2 * N)], temp)
                          
                          ###1e-6
                          start <- proc.time()[3]
                          Result <- vigor(Data$Y[1:N], ETA, Verbose = FALSE, RandomIni = TRUE, Thresholdvalue = 1e-6, Maxiteration = 10000)
                          end <- proc.time()[3]
                          
                          Result.sim2[35] <- end - start
                          Result.sim2[36] <- length(Result$ResidualVar)
                          
                          ####Estimation accuracy
                          Result.sim2[37] <- cor(Result$ETA[[1]]$Beta, as.numeric(Data$B1))
                          
                          ####AUC
                          ROCdata <- cbind(Result$ETA[[2]]$Rho, c(rep(1, p2*Nzero2), rep(0, p2-p2*Nzero2)))
                          ROCdata <- ROCdata[order(ROCdata[, 1], decreasing = TRUE), ]
                          Pred <- prediction(ROCdata[, 1], ROCdata[, 2])
                          AUC.temp <- performance(Pred, "auc")
                          Result.sim2[38] <- as.numeric(AUC.temp@y.values)
                          
                          ####Prediction accuracy
                          temp <- predict_vigor(Result, list(Data$X1[(N + 1):(2 * N), ], Data$X2[(N + 1):(2 * N), ]))
                          Result.sim2[39] <- cor(Data$U[(N + 1):(2 * N)], temp)
                          
                          Result.sim2
                        }#function
    )#parLapply
    
    temp <- matrix(unlist(Result),nr=39)
    Stats.Exp2.BRRBayesB <- rbind(Stats.Exp2.BRRBayesB, temp)
    
  }#block
}#i
stopCluster(cl)

#See the averages
for(i in 1:length(P1)){
  
  temp <- Stats.Exp2.BRRBayesB[((i - 1) * 39 * 4 + 1):(i * 39 * 4),]
  temp <- apply(matrix(as.numeric(temp), nr = 39), 1, mean)
  
  cat("#BGLR","\n")
  cat("#", P1[i] * 2, temp[1:8], "\n")
  cat("#", P1[i] * 2, temp[9:16], "\n")
  cat("#", P1[i] * 2, temp[17:24], "\n")
  cat("#VIGoR","\n")  
  cat("#", P1[i] * 2, temp[25:29], "\n")
  cat("#", P1[i] * 2, temp[30:34], "\n")
  cat("#", P1[i] * 2, temp[35:39], "\n\n")
}


#Exp. 3 (real data analysis)####################################################
##Real data is taken from Onogi et al. 2021 (https://www.frontiersin.org/articles/10.3389/fgene.2021.803636/full)
##Phenotype and meteorological information are provided in Supplementary Data 1
##which can be downloaded from https://doi.org/10.5061/dryad.rr4xgxd6r.
Data<-read.csv("SupplementaryData1.csv", header = TRUE)

#Use Variety V083 which has the most records
sort(table(Data$VarietyID))

Use.row<-Data$VarietyID == "V083" & !is.na(Data$DTF)
N <- sum(Use.row)
N
1051

#Use T (mean temperature), Ph (photo period length), and Pr (precipitation) as explanatory variables
Use.col <- c(grep("T[.]", colnames(Data)),
             grep("Ph[.]" , colnames(Data)),
             grep("Pr[.]" , colnames(Data)))
Use.col <- is.element(1:ncol(Data), Use.col)
sum(Use.col)
738
#=>246/meteorological factor x 3 = 738 explanatory variables

#Explanatory variable for additive effect
X0 <- Data[Use.row, Use.col]

#Impute explanatory variables with averages
for(i in 1:ncol(X0)){
  if(any(is.na(X0[, i])))
    X0[, i][is.na(X0[, i])] <- mean(X0[, i][!is.na(X0[, i])], na.rm = TRUE)
}

#Response variable (days to flowering)
Y <- Data$DTF[Use.row]
any(is.na(Y))
FALSE

#Standardize X0
X0 <- scale(X0)

#Explanatory variables for interaction effects
X1 <- X2 <- X3 <- matrix(0, N, 246 * 246)
#T x Ph
k <- 1
for(i in 1:246){
  for(j in 1:246){
    X1[, k] <- X0[, i] * X0[, j + 246]
    k <- k + 1
  }
}
#T x Pr
k <- 1
for(i in 1:246){
  for(j in 1:246){
    X2[, k] <- X0[, i] * X0[, j + 246 * 2]
    k <- k + 1
  }
}
#Ph x Pr
k <- 1
for(i in 1:246){
  for(j in 1:246){
    X3[, k] <- X0[, i + 246] * X0[, j + 246 * 2]
    k <- k + 1
  }
}

#Split data into training and test data
##data after 2005 is used for test
Train <- Data$Year[Use.row] < 2005
Test <- Data$Year[Use.row] >= 2005
sum(Train)
838
sum(Test)
213
rm(Data)

#Create hyperparameter values using hyperpara
H0 <- hyperpara(X0[Train, ], 1, "BayesC", 0.2, "Var")
H1 <- hyperpara(X1[Train, ], 1, "BayesC", 0.001, "Var")
H2 <- hyperpara(X2[Train, ], 1, "BayesC", 0.001, "Var")
H3 <- hyperpara(X3[Train, ], 1, "BayesC", 0.001, "Var")

#Run VIGoR
ETA <- list(list(model = "BayesC",
                 X = X0[Train, ],
                 H = H0),
            list(model = "BayesC",
                 X = X1[Train, ],
                 H = H1),
            list(model = "BayesC",
                 X = X2[Train, ],
                 H = H2),
            list(model = "BayesC",
                 X = X3[Train, ],
                 H = H3))

Predict.vigor.BayesC <- NULL
for(th in c(1e-4, 1e-5, 1e-6)){
  for(i in 1:10){
    cat(th, i,"\n")
    
    start <- proc.time()[3]
    Result.vigor <- vigor(Y[Train], ETA, Verbose = FALSE, Thresholdvalue = th, Maxiteration = 10000)
    end <- proc.time()[3]
    
    #Predict test data
    temp <- predict_vigor(Result.vigor, list(X0[Test, ], X1[Test, ], X2[Test, ], X3[Test, ]))
    
    Predict.vigor.BayesC <- rbind(Predict.vigor.BayesC,
                                  c(end - start, #calculation time
                                    length(Result.vigor$ResidualVar), #number of iterations
                                    cor(Y[Test], temp), #correlation
                                    sqrt(mean((Y[Test] - temp)^2)))) #RMSE
  }
}

##see the averages
for(i in 1:3){
  cat("#", c(1e-4, 1e-5, 1e-6)[i], "\n")
  cat("#", i, apply(Predict.vigor.BayesC[((i - 1) * 10 + 1):(i * 10), ], 2, mean),"\n")
  cat("#", i, apply(Predict.vigor.BayesC[((i - 1) * 10 + 1):(i *10), ], 2, sd),"\n")
}
rm(Result.vigor)


#Run BGLR
ETA <- list(list(model = "BayesC",
                 X = X0[Train, ],
                 df0 = H0[1], S0 = H0[2], probIn = H0[3], counts = 1e+6),
            list(model = "BayesC",
                 X = X1[Train, ],
                 df0 = H1[1], S0 = H1[2], probIn = H1[3], counts = 1e+6),
            list(model = "BayesC",
                 X = X2[Train, ],
                 df0 = H2[1], S0 = H2[2], probIn = H2[3], counts = 1e+6),
            list(model = "BayesC",
                 X = X3[Train, ],
                 df0 = H3[1], S0 = H3[2], probIn = H3[3], counts = 1e+6))

#MCMC conditions
##Three chain lengths
##Number of iterations, burnin, sampling interval, number of samples used for inference
Cond <- matrix(c(1500, 500, 5, 200,
                 15000, 5000, 10, 1000,
                 30000, 20000, 10, 1000),
               nc = 4, byrow = TRUE)

##Use snow
###10 replications are conducted using 5 clusters.
x<-1:10
cl <- makeCluster(5, type = "SOCK")
clusterExport(cl, c("BGLR", "Y", "Train", "Test", "ETA", "X0", "X1", "X2", "X3"))
Predict.BGLR.BayesC <- as.list(numeric(nrow(Cond)))

for(i in 1:nrow(Cond)){
  
  ite <- Cond[i, 1]
  burnin <- Cond[i, 2]
  th <- Cond[i, 3]
  sam <- Cond[i, 4]
  
  clusterExport(cl, c("ite", "burnin", "th", "sam"))
  cat(i, "\n")
  
  Result1 <- parSapply(cl,
                       x[1:5],
                       function(d){
                         
                         start <- proc.time()[3]
                         Result <- BGLR(y = Y[Train], ETA = ETA, verbose = FALSE,
                                        nIter = ite, burnIn = burnin, thin = th,
                                        saveAt = paste(d, "/", sep = ""))
                         end <- proc.time()[3]
                         
                         #Prediction accuracy
                         temp <- Result$mu + 
                           colSums(t(X0[Test, ]) * Result$ETA[[1]]$b) + 
                           colSums(t(X1[Test, ]) * Result$ETA[[2]]$b) + 
                           colSums(t(X2[Test, ]) * Result$ETA[[3]]$b) + 
                           colSums(t(X3[Test, ]) * Result$ETA[[4]]$b)
                         
                         #Rhat
                         s <- scan(paste(d, "varE.dat", sep = "/"))
                         s <- matrix(s[(length(s)-sam+1):length(s)], nc=2)
                         B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                         W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                         Rhat.varE <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                         
                         s <- read.table(paste(d, "ETA_1_parBayesC.dat", sep = "/"))[, 1]
                         s <- matrix(s[(length(s)-sam+1):length(s)], nc=2)
                         B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                         W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                         Rhat.varG1 <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                         
                         s <- read.table(paste(d, "ETA_2_parBayesC.dat", sep = "/"))[, 1]
                         s <- matrix(s[(length(s)-sam+1):length(s)], nc=2)
                         B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                         W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                         Rhat.varG2 <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                         
                         s <- read.table(paste(d, "ETA_3_parBayesC.dat", sep = "/"))[, 1]
                         s <- matrix(s[(length(s)-sam+1):length(s)], nc=2)
                         B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                         W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                         Rhat.varG3 <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                         
                         s <- read.table(paste(d, "ETA_4_parBayesC.dat", sep = "/"))[, 1]
                         s <- matrix(s[(length(s)-sam+1):length(s)], nc=2)
                         B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                         W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                         Rhat.varG4 <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                         
                         s <- scan(paste(d, "mu.dat", sep = "/"))
                         s <- matrix(s[(length(s)-sam+1):length(s)], nc=2)
                         B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                         W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                         Rhat.mu <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                         
                         #Stats
                         c(end - start,#calculation time
                           cor(Y[Test], temp),#Correlation
                           sqrt(mean((Y[Test] - temp)^2)),#RMSE
                           Rhat.varE,
                           Rhat.varG1,
                           Rhat.varG2,
                           Rhat.varG3,
                           Rhat.varG4,
                           Rhat.mu)
                       }
  )
  
  Result2 <- parSapply(cl,
                       x[6:10],
                       function(d){
                         
                         start <- proc.time()[3]
                         Result <- BGLR(y = Y[Train], ETA = ETA, verbose = FALSE,
                                        nIter = ite, burnIn = burnin, thin = th,
                                        saveAt = paste(d, "/", sep = ""))
                         end <- proc.time()[3]
                         
                         #Prediction accuracy
                         temp <- Result$mu + 
                           colSums(t(X0[Test, ]) * Result$ETA[[1]]$b) + 
                           colSums(t(X1[Test, ]) * Result$ETA[[2]]$b) + 
                           colSums(t(X2[Test, ]) * Result$ETA[[3]]$b) + 
                           colSums(t(X3[Test, ]) * Result$ETA[[4]]$b)
                         
                         #Rhat
                         s <- scan(paste(d, "varE.dat", sep = "/"))
                         s <- matrix(s[(length(s)-sam+1):length(s)], nc=2)
                         B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                         W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                         Rhat.varE <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                         
                         s <- read.table(paste(d, "ETA_1_parBayesC.dat", sep = "/"))[, 1]
                         s <- matrix(s[(length(s)-sam+1):length(s)], nc=2)
                         B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                         W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                         Rhat.varG1 <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                         
                         s <- read.table(paste(d, "ETA_2_parBayesC.dat", sep = "/"))[, 1]
                         s <- matrix(s[(length(s)-sam+1):length(s)], nc=2)
                         B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                         W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                         Rhat.varG2 <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                         
                         s <- read.table(paste(d, "ETA_3_parBayesC.dat", sep = "/"))[, 1]
                         s <- matrix(s[(length(s)-sam+1):length(s)], nc=2)
                         B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                         W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                         Rhat.varG3 <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                         
                         s <- read.table(paste(d, "ETA_4_parBayesC.dat", sep = "/"))[, 1]
                         s <- matrix(s[(length(s)-sam+1):length(s)], nc=2)
                         B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                         W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                         Rhat.varG4 <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                         
                         s <- scan(paste(d, "mu.dat", sep = "/"))
                         s <- matrix(s[(length(s)-sam+1):length(s)], nc=2)
                         B <- sum((apply(s, 2, mean) - mean(s))^2) * nrow(s) / (2 - 1)
                         W <- sum((t(s) - apply(s, 2, mean))^2)/((nrow(s) - 1) * 2)
                         Rhat.mu <- sqrt((W * (nrow(s) - 1)/nrow(s) + B / nrow(s)) / W)
                         
                         #Stats
                         c(end - start,
                           cor(Y[Test], temp),
                           sqrt(mean((Y[Test] - temp)^2)),
                           Rhat.varE,
                           Rhat.varG1,
                           Rhat.varG2,
                           Rhat.varG3,
                           Rhat.varG4,
                           Rhat.mu)
                       }
  )
  
  Predict.BGLR.BayesC[[i]] <-cbind(Result1, Result2)
}
stopCluster(cl)

##See the averages
for(i in 1:nrow(Cond)){
  cat("#", Cond[i, 1], "\n")
  cat("#", i, apply(Predict.BGLR.BayesC[[i]], 1, mean), "\n")
  cat("#", i, apply(Predict.BGLR.BayesC[[i]], 1, sd), "\n")
}

