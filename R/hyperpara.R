hyperpara <- function(X, Mvar,
                      Model = c("BL", "EBL", "BayesA", "BayesB", "BayesC", "BRR"),
                      Kappa = 0.01, Xtype = c("Geno", "Var"), f = 0, BL.Phi = 1, EBL.Phi = 0.1,
                      EBL.Omega = 0.1, EBL.Psi = 1, Nu = 5, Verbose=FALSE){

  stopifnot(is.matrix(X))
  if(any(is.na(X))) stop("NA in X is not allowed")
  P <- ncol(X)
  N <- nrow(X)

  Model <- match.arg(Model)
  if(Model == "BL" | Model == "EBL"){
    if(any(Mvar >=1 | Mvar <= 0))
      stop("Mvar should be 0<Mvar<1 when BL or EBL")
  } else {
    if(any(Mvar > 1 | Mvar <= 0))
      stop("Mvar should be 0<Mvar<=1")
  }

  if(any(Kappa > 1 | Kappa <= 0))
    stop("Kappa should be 0<Kappa<=1")
  if(Xtype[1] != "Geno" & Xtype[1] != "Var")
    stop("Xtype specification error")
  if((Xtype[1] == "Geno") & any(X > 2 | X < 0))
    stop("Genotypes should be coded between 0 and 2")
  if(length(f) > 1 | f < 0 | f > 1)
    stop("f should be a scalar (0 <= f <= 1)")
  if(any(BL.Phi <= 0))
    stop("BL.Phi should be >0")
  if(any(EBL.Phi <= 0))
    stop("EBL.Phi should be >0")
  if(any(EBL.Omega <= 0))
    stop("EBL.Omega should be >0")
  if(any(EBL.Psi < 0))
    stop("EBL.Psi should be >=0")
  if(any(Nu <= 2))
    stop("Nu should be >2")

  if(Xtype[1] == "Var"){
    Sum2pq <- sum(apply(X, 2, var))
    if(Verbose){
      cat("\n")
      cat("Model:", Model, "\n")
      cat("Mvar:", Mvar, "\n")
      cat("N. of individuals:", N, "\n")
      cat("N. of markers:", P, "\n")
    }
  }else{
    Af <- colSums(X)/2/N
    Sum2pq <- sum(2 * (1 + f) * Af * (1 - Af))
    Af[Af > 0.5] <- 1 - Af[Af > 0.5]
    if(Verbose){
      cat("\n")
      cat("Model:", Model, "\n")
      cat("Mvar:", Mvar, "\n")
      cat("N. of individuals:", N, "\n")
      cat("N. of markers:", P, "\n")
    }
  }

  if(Model == "BL"){
    L.Mvar <- length(Mvar)
    L.Kappa <- length(Kappa)
    L.Phi <- length(BL.Phi)
    Nset <- L.Mvar * L.Kappa * L.Phi
    Comb <-cbind(rep(Mvar, each = L.Kappa * L.Phi),
                 rep(Kappa, each = L.Phi),
                 rep(BL.Phi))

    HM <- cbind(Comb[, 3],
                Comb[, 3]/(2 * Comb[, 2] * Sum2pq * (1/Comb[, 1] - 1)))

    if(Nset == 1){
      HM <- as.vector(HM); names(HM) <- c("Phi", "Omega")
    } else {
      colnames(HM) <- c("Phi", "Omega")
      rownames(HM) <- paste("Mvar", Comb[, 1],
                            "_Kappa", Comb[, 2],
                            "_Phi", Comb[, 3], sep="")

    }
  }
  if(Model == "EBL"){
    L.Mvar <- length(Mvar)
    L.Kappa <- length(Kappa)
    L.Phi <- length(EBL.Phi)
    L.Omega <- length(EBL.Omega)
    L.EBL.Psi <- length(EBL.Psi)
    Nset <- L.Mvar * L.Kappa * L.Phi * L.Omega * L.EBL.Psi
    Comb <- cbind(rep(Mvar, each = L.Kappa * L.Phi * L.Omega * L.EBL.Psi),
                  rep(Kappa, each = L.Phi * L.Omega * L.EBL.Psi),
                  rep(EBL.Phi, each = L.Omega * L.EBL.Psi),
                  rep(EBL.Omega, each = L.EBL.Psi),
                  rep(EBL.Psi))

    HM <- cbind(Comb[, 3],
                Comb[ ,4],
                Comb[ ,5],
                Comb[, 3]/Comb[, 4] *
                  Comb[, 5]/(2 * Comb[, 2] * Sum2pq * (1/Comb[, 1] - 1)))

    if(Nset == 1){
      HM <- as.vector(HM); names(HM) <- c("Phi", "Omega", "Psi", "Theta")
    } else {
      colnames(HM)<-c("Phi","Omega","Psi","Theta")
      rownames(HM) <- paste("Mvar", Comb[, 1],
                            "_Kappa", Comb[, 2],
                            "_Phi", Comb[, 3],
                            "_Omega", Comb[, 4],
                            "_Psi", Comb[, 5], sep="")
    }
  }
  if(Model == "BayesA" | Model == "BRR"){
    Kappa <- 1
    L.Mvar <- length(Mvar)
    L.Nu <- length(Nu)
    Nset <- L.Mvar * L.Nu
    Comb <-cbind(rep(Mvar, each = L.Nu),
                 rep(Nu))
    HM <- cbind(Comb[, 2],
                (Comb[, 2] - 2)/Comb[, 2] * Comb[, 1]/Sum2pq)

    if(Nset == 1){
      HM <- as.vector(HM); names(HM) <- c("Nu", "S2")
    }else{
      colnames(HM) <- c("Nu", "S2")
      rownames(HM) <- paste("Mvar", Comb[, 1],
                            "_Nu", Comb[, 2], sep="")
    }
  }
  if(Model == "BayesC" | Model == "BayesB"){
    L.Mvar <- length(Mvar)
    L.Kappa <- length(Kappa)
    L.Nu <- length(Nu)
    Nset <- L.Mvar * L.Kappa * L.Nu
    Comb <-cbind(rep(Mvar, each = L.Kappa * L.Nu),
                 rep(Kappa, each = L.Nu),
                 rep(Nu))
    HM <- cbind(Comb[, 3],
                (Comb[, 3] - 2)/Comb[, 3] * Comb[, 1]/(Comb[, 2] * Sum2pq),
                Comb[, 2])

    if(Nset == 1){
      HM <- as.vector(HM); names(HM) <- c("Nu", "S2", "Kappa")
    }else{
      colnames(HM) <- c("Nu", "S2", "Kappa")
      rownames(HM) <- paste("Mvar", Comb[, 1],
                            "_Kappa", Comb[, 2],
                            "_Nu", Comb[, 3], sep="")
    }
  }
  HM
}
