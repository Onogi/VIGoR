vigor <- function(Y, ETA, Function = c("fitting", "tuning", "cv"), Nfold = 5,
                  CVFoldTuning = 5, Partition = NULL, Thresholdvalue = 1e-5,
                  Maxiteration = 1000, RandomIni = TRUE,
                  Metrics = c("rmse", "cor"), Verbose = TRUE){

  #Functions
  genomewideregression <- function(y, eta, priortype, methodcode, p, h,
                                   thresholdvalue, maxiteration, randomize,
                                   variance, division_H, division_V){

    #y and x should be sorted in the same order
    nm <- length(priortype)
    n <- length(y)
    ytox <- 0:(n - 1)
    xtoy <- 0:(n - 1)
    tau0 <- c(1, 0)
    lb <- 0
    #lb <- numeric(maxiteration)
    rmonitor <- numeric(maxiteration)

    #standardize y
    y.mean <- mean(y)
    y.sd <- sd(y)
    y.scaled <- (y - y.mean)/y.sd

    #Concatenate X
    x <- NULL
    for(method in 1:nm)
      if(eta[[method]]$model == "BLUP"){
        x <- cbind(x,
                   eta[[method]]$values,
                   eta[[method]]$vectors,
                   eta[[method]]$tvectors)
      }else{
        x <- cbind(x, eta[[method]]$X)
      }

    #Objects for storing results
    expectation <- NULL
    uncertainty <- NULL
    gamma <- NULL
    for(method in 1:nm)
      if(eta[[method]]$model == "BLUP") {
        expectation <- c(expectation, numeric(n))
        uncertainty <- c(uncertainty, numeric(n))
        gamma <- c(gamma, numeric(n))
      }else {
        expectation <- c(expectation, numeric(p[method]))
        uncertainty <- c(uncertainty, numeric(p[method]))
        gamma <- c(gamma, numeric(p[method]))
      }

    #Divisions of concatenated objects
    division_G <- numeric(nm)
    division_E <- numeric(nm)
    if(nm > 1)
      for(method in 2:nm)
        if(eta[[method-1]]$model == "BLUP"){
          division_G[method] <- division_G[method-1] + 1 + 2 * n
          division_E[method] <- division_E[method-1] + n
        }else{
          division_G[method] <- division_G[method-1] + p[method-1]
          division_E[method] <- division_E[method-1] + p[method-1]
        }
    division_G <- n * division_G

    #Fitted values
    fittedvalue <- numeric(n)

    #Run
    Result<-.C("vigor_c",
               as.integer(priortype),
               as.integer(methodcode),
               as.integer(nm),
               as.integer(p),
               as.integer(n),
               as.integer(n),
               as.integer(ytox),
               as.integer(xtoy),
               as.integer(maxiteration),
               as.integer(randomize),
               as.integer(division_G),
               as.integer(division_H),
               as.integer(division_E),
               as.integer(division_V),
               as.double(thresholdvalue),
               as.double(y.scaled),
               as.double(x),#17
               as.double(h),#18
               as.double(tau0),#19
               as.double(lb),#20
               as.double(rmonitor),#21
               as.double(expectation),#22
               as.double(uncertainty),#23
               as.double(variance),#24
               as.double(gamma),
               as.double(fittedvalue),
               PACKAGE="VIGoR"
    )

    #Extract results
    Result_output <- list(LB = Result[[20]],
                          ResidualVar = Result[[21]] * y.sd^2,
                          H = h,
                          Fittedvalue = Result[[26]] * y.sd + y.mean)
    Result_output$ResidualVar <-
      Result_output$ResidualVar[Result_output$ResidualVar > 0]

    Result_ETA <- as.list(numeric(nm))
    for(method in 1:nm){
      s1 <- division_E[method] + 1
      if(method == nm)
        e1 <- length(Result[[22]])
      else
        e1 <- division_E[method + 1]

      s2 <- division_V[method] + 1
      if(method == nm)
        e2 <- length(Result[[24]])
      else
        e2 <- division_V[method + 1]

      #FIXED
      if(methodcode[method] == 9){
        Result_ETA[[method]] <- list(Beta    = Result[[22]][s1:e1] * y.sd,
                                     Sd.beta = Result[[23]][s1:e1] * y.sd)
        #The first column is regarded as the intercept
        Result_ETA[[method]]$Beta[1] <- Result_ETA[[method]]$Beta[1] + y.mean
      }

      #BL and EBL
      if(methodcode[method] == 1|methodcode[method] == 2)
        Result_ETA[[method]] <- list(Beta    = Result[[22]][s1:e1] * y.sd,
                                     Sd.beta = Result[[23]][s1:e1] * y.sd)
      #BayesA, BayesB, BayesC, and BRR
      if(methodcode[method] == 4|methodcode[method] == 7)
        if(h[division_H[method] + 3] == 1){
          Result_ETA[[method]] <- list(Beta    = Result[[22]][s1:e1] * y.sd,
                                       Sd.beta = Result[[23]][s1:e1] * y.sd,
                                       Sigma2  = Result[[24]][s2:e2] * y.sd^2)
        }else{
          Result_ETA[[method]] <- list(Beta    = Result[[22]][s1:e1] * y.sd,
                                       Sd.beta = Result[[23]][s1:e1] * y.sd,
                                       Sigma2  = Result[[24]][s2:e2] * y.sd^2,
                                       Rho    = Result[[25]][s1:e1])
        }
      #BLUP
      if(methodcode[method] == 8)
        Result_ETA[[method]] <- list(U      = Result[[22]][s1:e1] * y.sd,
                                     Sd.u   = Result[[23]][s1:e1] * y.sd,
                                     Sigma2 = Result[[24]][s2:e2] * y.sd^2,
                                     iK = eta[[method]]$iK)
    }
    Result_output <- c(Result_output, ETA = list(Result_ETA))
    Result_output
  }

  tuning <- function(y, eta, priortype, methodcode, p, h,
                     thresholdvalue, maxiteration, randomize, cvfoldtuning,
                     variance, division_H, division_V, objfunction, verbose){

    n <- length(y)
    nset <- nrow(h)
    obj <- numeric(nset)
    nm <- length(methodcode)

    am <- n %% cvfoldtuning
    if(am == 0) n2 <- n else n2 <- n + cvfoldtuning - am
    shuffled <- matrix(sample(1:n2, n2, replace=FALSE), ncol = cvfoldtuning)
    test <- as.list(numeric(cvfoldtuning))
    for(fold in 1:cvfoldtuning)
      test[[fold]] <- sort(shuffled[, fold][shuffled[, fold] <= n])

    for(set in 1:nset){
      prediction <- numeric(n)
      for(fold in 1:cvfoldtuning){
        if(verbose) cat("tuning", "set", set, "fold", fold, "\n")
        eta.fold <- eta
        for(method in 1:nm){
         if(eta.fold[[method]]$model == "BLUP"){
           K.eigen <-
             eigen(eta[[method]]$K[-test[[fold]], -test[[fold]], drop = FALSE])
           iK <- K.eigen$vectors %*%
             diag(1/K.eigen$values) %*%
             t(K.eigen$vectors)
           eta.fold[[method]]$iK <- iK
           eta.fold[[method]]$values <- 1/K.eigen$values
           eta.fold[[method]]$vectors <- K.eigen$vectors
           eta.fold[[method]]$tvectors <- t(K.eigen$vectors)
         }else{
           eta.fold[[method]]$X <-
             eta.fold[[method]]$X[-test[[fold]], , drop = FALSE]
         }
        }
        result <- genomewideregression(y[-test[[fold]]],
                                       eta.fold,
                                       priortype,
                                       methodcode,
                                       p,
                                       h[set,],
                                       thresholdvalue,
                                       maxiteration,
                                       randomize,
                                       variance,
                                       division_H,
                                       division_V)
        for(method in 1:nm){
          if(methodcode[method] == 8){
            #BLUP
            prediction[test[[fold]]] <-
              prediction[test[[fold]]] +
              eta[[method]]$K[test[[fold]], -test[[fold]]] %*%
              eta.fold[[method]]$iK %*%
              matrix(result$ETA[[method]]$U, ncol = 1)
          }else{
            prediction[test[[fold]]] <-
              prediction[test[[fold]]] +
              eta[[method]]$X[test[[fold]], , drop=FALSE] %*%
              matrix(result$ETA[[method]]$Beta, ncol = 1)
          }
        }
      }
      obj[set] <- switch(objfunction,
                         rmse = sqrt(sum((prediction - y)^2, na.rm=T)/N),
                         cor = cor(as.numeric(prediction), as.numeric(y)))
    }
    obj
  }

  #Add intercept when no fixed effect is included
  stopifnot(is.vector(Y))
  Nall <- length(Y)
  Nm <- length(ETA)
  Nfixed <- 0
  for(method in 1:Nm)
    if(ETA[[method]]$model == "FIXED") Nfixed <- Nfixed + 1
  if(Nfixed > 1)
    stop("Fixed effects should be given in a single method")
  AddIntercept <- FALSE
  if(Nfixed == 0) AddIntercept <- TRUE
  if(AddIntercept){
    ETA <- c(ETA,
             list(list(model = "FIXED", X = matrix(1, nrow = Nall, ncol = 1))))
    Nm <- Nm + 1
  }

  #Check ETA
  for(method in 1:Nm)
    if(is.null(ETA[[method]]$model))
      stop("Please specify model for method ", method)

  for(method in 1:Nm)
    if(is.element(ETA[[method]]$model,
                  c("BL", "EBL", "BayesA", "BayesB", "BayesC", "BRR"))){
      #X should be included
      if(!any(names(ETA[[method]]) == "X"))
        stop("X should be included in ETA for method ", method)
      if(!is.matrix(ETA[[method]]$X))
        stop("X should be a matrix (method ", method, ")")
      if(nrow(ETA[[method]]$X) != Nall)
        stop("The number of rows of X is inconsistent with Y (method ",
             method, ")")

    }else{
      if(ETA[[method]]$model == "BLUP"){
        #Either X or K should be included
        if(!any(is.element(names(ETA[[method]]), c("X", "K"))))
          stop("Either X or K should be included in ETA for method ", method)
        if(any(names(ETA[[method]]) == "X")) {
          if(!is.matrix(ETA[[method]]$X))
            stop("X should be a matrix (method ", method, ")")
          if(nrow(ETA[[method]]$X) != Nall)
            stop("The number of rows of X is inconsistent with Y (method ",
                 method, ")")
        }
        if(any(names(ETA[[method]]) == "K")) {
          if(!is.matrix(ETA[[method]]$K))
            stop("K should be a matrix (method ", method, ")")
          if(nrow(ETA[[method]]$K) != Nall | ncol(ETA[[method]]$K) != Nall)
            stop("The dimensions of K is inconsistent with Y (method ",
                 method, ")")
        }

      }else{
        if(ETA[[method]]$model == "FIXED"){
          #Two kinds of specification
          #Using formula with data or using matrix (X)
          Class <- lapply(ETA[[method]], class)
          if(any(Class == "formula")) {
            if(!any(names(ETA[[method]]) == "data"))
              stop("data is not included in ETA for method ", method,
                   " (only formula found)")
            if(nrow(ETA[[method]]$data) != Nall)
              stop("The number of rows of data is inconsistent with Y (method ",
                   method, ")")
          }else{
            if(!any(names(ETA[[method]]) == "X"))
              stop("X should be included in ETA for method ", method)
            if(!is.matrix(ETA[[method]]$X))
              stop("X should be a matrix (method ", method, ")")
            if(nrow(ETA[[method]]$X) != Nall)
              stop("The number of rows of X is inconsistent with Y (method ",
                   method, ")")
          }
        }else{
          stop("Misspecification of model for method ", method)
        }
      }
    }

  #Check the number of samples
  for(method in 1:Nm)
    if(is.element(ETA[[method]]$model,
                  c("BL", "EBL", "BayesA", "BayesB", "BayesC", "BRR"))){
      if(nrow(ETA[[method]]$X) != Nall)
        stop("nrow(X) for method ", method, " is inconsistent with length(Y)")
    }else{
      if(ETA[[method]]$model == "BLUP"){
        if(any(names(ETA[[method]]) == "K")){
          if(ncol(ETA[[method]]$K) != nrow(ETA[[method]]$K))
            stop("Row and column numbers of K are inconsistent (method ",
                 method, ")")
          if(nrow(ETA[[method]]$K) != Nall)
            stop("nrow(K) for method ", method,
                 " is inconsistent with length(Y)")
        }else{
          if(nrow(ETA[[method]]$X) != Nall)
            stop("nrow(X) for method ", method,
                 " is inconsistent with length(Y)")
        }
      }else{
        #FIXED
        Class <- lapply(ETA[[method]], class)
        if(any(Class == "formula")) {
          if(nrow(ETA[[method]]$data) != Nall)
            stop("nrow(data) for method ", method,
                 " is inconsistent with length(Y)")
        }else{
          if(nrow(ETA[[method]]$X) != Nall)
            stop("nrow(X) for method ", method,
                 " is inconsistent with length(Y)")
        }
      }
    }

  #Use samples with response variables (Y)
  Use <- !is.na(Y)
  N <- sum(Use)
  Y <- Y[Use]
  for(method in 1:Nm)
    if(is.element(ETA[[method]]$model,
                  c("BL", "EBL", "BayesA", "BayesB", "BayesC", "BRR"))){
      ETA[[method]]$X <- ETA[[method]]$X[Use, , drop = FALSE]
      if(any(is.na(ETA[[method]]$X)))
        stop("NA in X is not allowed (method ", method, ")")
    }else{
      if(ETA[[method]]$model == "BLUP"){
        if(any(names(ETA[[method]]) == "K")){
          ETA[[method]]$K <- ETA[[method]]$K[Use, Use, drop = FALSE]
          if(any(is.na(ETA[[method]]$K)))
            stop("NA in K is not allowed (method ", method, ")")
        }else{
          ETA[[method]]$X <- ETA[[method]]$X[Use, , drop = FALSE]
          if(any(is.na(ETA[[method]]$X)))
            stop("NA in X is not allowed (method ", method, ")")
        }
      }else{
        #FIXED
        Class <- lapply(ETA[[method]], class)
        if(any(Class == "formula")) {
          ETA[[method]]$data <- ETA[[method]]$data[Use, ,drop = FALSE]
          if(any(is.na(ETA[[method]]$data)))
            stop("NA in data is not allowed (method ", method, ")")
        }else{
          ETA[[method]]$X <- ETA[[method]]$X[Use, , drop = FALSE]
          if(any(is.na(ETA[[method]]$X)))
            stop("NA in X is not allowed (method ", method, ")")
        }
      }
    }

  #Check K for BLUP
  CreateLinearKernel <- rep(FALSE, Nm)
  for(method in 1:Nm)
    if(ETA[[method]]$model == "BLUP"){
      if(!any(names(ETA[[method]]) == "K")){
        #Linear kernel is created from X
        X <- scale(ETA[[method]]$X)
        K <- X %*% t(X) / ncol(X)
        ETA[[method]] <- c(ETA[[method]], K = list(K))
        CreateLinearKernel[method] <- TRUE
        rm(X, K)
      }
      K.eigen <- eigen(ETA[[method]]$K)
      if(any(K.eigen$values <= 0)) {
        diag(ETA[[method]]$K) <- diag(ETA[[method]]$K) +
          mean(diag(ETA[[method]]$K)) * 1e-3
        K.eigen <- eigen(ETA[[method]]$K)
        if(any(K.eigen$values <= 0))
          stop ("Eigen values of K include values <= 0 (method ", method, ")")
      }
      ETA[[method]] <- c(ETA[[method]],
                         iK = list(K.eigen$vectors %*%
                                     diag(1/K.eigen$values) %*%
                                     t(K.eigen$vectors)),
                         values = list(1/K.eigen$values),
                         vectors = list(K.eigen$vectors),
                         tvectors = list(t(K.eigen$vectors)))
      rm(K.eigen)
    }

  #Transform factors to model matrices
  for(method in 1:Nm)
    if(ETA[[method]]$model == "FIXED"){
      Class <- lapply(ETA[[method]], class)
      if(any(Class == "formula")) {
        Which <- which(Class == "formula")[1]
        X <- model.matrix(ETA[[method]][[Which]], ETA[[method]]$data)
        if(!is.null(ETA[[method]]$X)){
          #overwrite
          ETA[[method]]$X <- X
        }else{
          ETA[[method]] <- c(ETA[[method]], X = list(X))
        }
      }
    }

  #Check hyperparameters
  Nh <- c(2, 4, 2, 3, 3, 2, 2, 0)#number of hyperparameters for each method
  names(Nh) <- c("BL", "EBL", "BayesA", "BayesB", "BayesC",
                 "BRR", "BLUP", "FIXED")
  Nset <- numeric(Nm)

  for(method in 1:Nm){
    if(is.null(ETA[[method]]$H)){
      #Assign default values
      ETA[[method]]$H <- switch(ETA[[method]]$model,
                                BL     = matrix(c(1, 1), nrow=1),
                                EBL    = matrix(c(0.1, 0.1, 1, 0.1), nrow=1),
                                BayesA = matrix(c(5, 0.01, 1), nrow=1),
                                BayesB = matrix(c(5, 0.1, 0.01), nrow=1),
                                BayesC = matrix(c(5, 0.1, 0.01), nrow=1),
                                BRR    = matrix(c(5, 0.01, 1), nrow=1),
                                BLUP   = matrix(c(5, 0.3), nrow=1))
    }else{
      #Check given values
      if(is.matrix(ETA[[method]]$H)){
        if(ncol(ETA[[method]]$H) != Nh[ETA[[method]]$model])
          stop("Number of hyperperameters is incorrect (method ", method, ")")
      }else{
        H <- as.vector(ETA[[method]]$H)
        if(length(H) != Nh[ETA[[method]]$model])
          stop("Number of hyperperameters is incorrect (method ", method, ")")
        ETA[[method]]$H <- matrix(H, nrow = 1)
      }

      if(ETA[[method]]$model == "BL"|ETA[[method]]$model == "EBL")
        if(any(ETA[[method]]$H <= 0))
          stop("Hyperparameters should be positive (method ", method, ")")
      if(ETA[[method]]$model == "BayesA"){
        if(any(ETA[[method]]$H[,1] <= 0))
          stop("Nu should be >0 (method ", method, ")")
        if(any(ETA[[method]]$H[,2] < 0))
          stop("S2 should be >=0 (method ", method, ")")
        #Add kappa
        ETA[[method]]$H <- cbind(ETA[[method]]$H, 1)
        colnames(ETA[[method]]$H)[3] <- "Kappa"
      }
      if(ETA[[method]]$model == "BayesB"){
        if(any(ETA[[method]]$H[,1] <= 0))
          stop("Nu should be >0 (method ", method, ")")
        if(any(ETA[[method]]$H[,2] < 0))
          stop("S2 should be >=0 (method ", method, ")")
        if(any(ETA[[method]]$H[,3] > 1|ETA[[method]]$H[,3] <= 0))
          stop("Kappa should be 0<Kappa<=1 (method ", method, ")")
      }
      if(ETA[[method]]$model == "BayesC"){
        if(any(ETA[[method]]$H[,1] <= 0))
          stop("Nu should be >0 (method ", method, ")")
        if(any(ETA[[method]]$H[,2] < 0))
          stop("S2 should be >=0 (method ", method, ")")
        if(any(ETA[[method]]$H[,3] > 1|ETA[[method]]$H[,3] <= 0))
          stop("Kappa should be 0<Kappa<=1 (method ", method, ")")
      }
      if(ETA[[method]]$model == "BRR"){
        if(any(ETA[[method]]$H[,1] <= 0))
          stop("Nu should be >0 (method ", method, ")")
        if(any(ETA[[method]]$H[,2] < 0))
          stop("S2 should be >=0 (method ", method, ")")
        #Add kappa
        ETA[[method]]$H <- cbind(ETA[[method]]$H, 1)
        colnames(ETA[[method]]$H)[3] <- "Kappa"
      }
      if(ETA[[method]]$model == "BLUP"){
        if(any(ETA[[method]]$H[,2] < 0))
          stop("S2 should be >=0 (method ", method, ")")
      }
      Nset[method] <- nrow(ETA[[method]]$H)
    }
    #Name hyperparameters
    colnames(ETA[[method]]$H) <-
      switch(ETA[[method]]$model,
             BL     = c("Phi", "Omega"),
             EBL    = c("Phi", "Omega", "Psi", "Theta"),
             BayesA = c("Nu", "S2", "Kappa"),
             BayesB = c("Nu", "S2", "Kappa"),
             BayesC = c("Nu", "S2", "Kappa"),
             BRR    = c("Nu", "S2", "Kappa"),
             BLUP   = c("Nu", "S2"))
  }

  #Concatenate hyperparameters of all methods
  Hyperparameters <- NULL
  Names <- NULL
  Where <- which(Nset > 0)
  Nsetforprod <- Nset
  Nsetforprod[Nsetforprod == 0] <- 1

  if(any(Nset > 1)){
    Nset.all <- prod(Nsetforprod)
    if(length(Where) > 1)
      for(method in Where[-length(Where)]){
        Each <- Nset.all/prod(Nsetforprod[Where[1]:Where[Where==method]])
        temp <- ETA[[method]]$H[rep(1:Nset[method], each=Each), ]
        for(i in 1:ncol(temp)){
          Hyperparameters <- cbind(Hyperparameters, temp[,i])
          Names <- c(Names, colnames(temp)[i])
        }
      }
    temp <- ETA[[Where[length(Where)]]]$H
    for(i in 1:ncol(temp)){
      Hyperparameters <- cbind(Hyperparameters, temp[,i])
      Names <- c(Names, colnames(temp)[i])
    }
  }else{
    Nset.all <- 1
    for(method in 1:Nm){
      Hyperparameters <- cbind(Hyperparameters, ETA[[method]]$H)
      Names <- c(Names, colnames(ETA[[method]]$H))
    }
  }
  colnames(Hyperparameters) <- Names

  #Create objects used in C source codes
  Priortype <- numeric(Nm)
  Methodcode <- numeric(Nm)
  for(method in 1:Nm){
    Priortype[method] <- switch(ETA[[method]]$model,
                                BL     = 1,
                                EBL    = 1,
                                BayesA = 2,
                                BayesB = 2,
                                BayesC = 2,
                                BRR    = 2,
                                BLUP   = 3,
                                FIXED  = 4)
    Methodcode[method] <- switch(ETA[[method]]$model,
                                 BL     = 1,
                                 EBL    = 2,
                                 BayesA = 7,
                                 BayesB = 7,
                                 BayesC = 4,
                                 BRR    = 4,
                                 BLUP   = 8,
                                 FIXED  = 9)
  }

  P <- numeric(Nm)
  for(method in 1:Nm)
    if(ETA[[method]]$model == "BLUP") {
      P[method] <- 0
    } else {
      P[method] <- ncol(ETA[[method]]$X)
    }

  Variance <- NULL
  for(method in 1:Nm){
    if(ETA[[method]]$model == "BayesA" | ETA[[method]]$model == "BayesB")
      Variance <- c(Variance, numeric(P[method]))
    if(ETA[[method]]$model == "BayesC" |
       ETA[[method]]$model == "BRR" |
       ETA[[method]]$model == "BLUP")
      Variance <- c(Variance, 0)
  }

  #Divisions of concatenated objects
  Division_H <- numeric(Nm)
  Division_V <- numeric(Nm)
  if(Nm > 1)
    for(method in 2:Nm){
      temp <- ncol(ETA[[method-1]]$H)
      if(!is.null(temp)){
        Division_H[method] <- Division_H[method-1] + temp
      }else{
        Division_H[method] <- Division_H[method-1]
      }
      temp <- switch(ETA[[method-1]]$model,
                     BL     = 0,
                     EBL    = 0,
                     BayesA = P[method-1],
                     BayesB = P[method-1],
                     BayesC = 1,
                     BRR    = 1,
                     BLUP   = 1,
                     FIXED  = 0)
      Division_V[method] <- Division_V[method-1] + temp
    }

  #Check Function
  if(!is.element(Function[1], c("fitting", "tuning", "cv")))
    stop("Misspecification of Function")
  if(Function[1] == "cv"){
    if(is.null(Partition)){
      if(Nfold > 1){
        Nfold <- round(Nfold)
        Am <- N %% Nfold
        if(Am == 0) N2 <- N else N2 <- N + Nfold - Am
        Partition <- matrix(sample(1:N2, N2, replace = FALSE), ncol = Nfold)
        Partition[Partition > N] <- -9
        Random <- TRUE
      }
      if(Nfold == -1){
        Partition<-matrix(1:N, ncol = N)
        Nfold <- N
        Random <- FALSE
      }
      if((Nfold <= 1&Nfold != -1)|Nfold > N)
        stop ("Nfold specification error")
    }else{
      if(!is.matrix(Partition)|!is.numeric(Partition) |
         any(is.na(Partition))|any(Partition == 0, na.rm = TRUE))
        stop("Partition matrix error")
      if(any(Partition > N, na.rm = TRUE))
        stop("Some numbers in Partition are larger than the length of Y")
      Nfold <- ncol(Partition)
      Random <- FALSE
    }
  }else{
    if(Function[1] != "fitting"&Function[1] != "tuning")
      stop("Function specification error")
  }

  #Check the maximum number of iterations
  stopifnot(Maxiteration > 0)

  #Set Randomize
  if(RandomIni) {Randomize <- 1} else {Randomize <- 0}

  #Check objective function of cross-validation
  if(!is.element(Metrics[1], c("rmse", "cor")))
    stop("Misspecification of Metrics")

  #Print information
  if(Verbose){
    cat("\n")
    cat("#Number of used samples:",N,"\n")
    cat("#Regression methods and hyperparameters:\n")
    for(method in 1:Nm){
      cat("Method", method, ":", ETA[[method]]$model, "\n")
      cat(switch(ETA[[method]]$model,
                 BL     = paste("Phi","Omega\n"),
                 EBL    = paste("Phi","Omega","Psi","Theta\n"),
                 BayesA = paste("v", "S2", "Kappa\n"),
                 BayesB = paste("v", "S2", "Kappa\n"),
                 BayesC = paste("v", "S2", "Kappa\n"),
                 BRR    = paste("v", "S2", "Kappa\n"),
                 BLUP   = paste("v", "S2\n")
                 ))
      if(ETA[[method]]$model != "FIXED")
        for(set in 1:Nset[method])
          cat(ETA[[method]]$H[set, ],"\n")
      if(ETA[[method]]$model != "BLUP")
        cat("Number of covariates:", P[method], "\n")
      if(method == Nm&AddIntercept)
        cat("(Intercept automatically added)\n")
      if(CreateLinearKernel[method])
        cat("Linear kernel was created from X\n")
      cat("\n")
    }
  }

  #Calculation
  if(Function[1] == "fitting"){
    if(Verbose) {
      cat("Model fitting\n")
      if(Nset.all > 1)
        cat("Multiple hyperparameter sets are found. The first one is used\n")
    }
    Result <- genomewideregression(Y,
                                   ETA,
                                   Priortype,
                                   Methodcode,
                                   P,
                                   Hyperparameters[1,],
                                   Thresholdvalue,
                                   Maxiteration,
                                   Randomize,
                                   Variance,
                                   Division_H,
                                   Division_V)
  }

  if(Function[1] == "tuning"){
    if(Nset.all == 1){
      if(Verbose)
        cat("Multiple hyperparameter sets are not found. Conduct fitting\n")
      Result <- genomewideregression(Y,
                                     ETA,
                                     Priortype,
                                     Methodcode,
                                     P,
                                     Hyperparameters,
                                     Thresholdvalue,
                                     Maxiteration,
                                     Randomize,
                                     Variance,
                                     Division_H,
                                     Division_V)
    }else{
      if(Verbose)
        cat("Model fitting after hyperparameter tuning\n")
      Obj <- tuning(Y,
                    ETA,
                    Priortype,
                    Methodcode,
                    P,
                    Hyperparameters,
                    Thresholdvalue,
                    Maxiteration,
                    Randomize,
                    CVFoldTuning,
                    Variance,
                    Division_H,
                    Division_V,
                    Metrics[1],
                    Verbose)
      Bestset <- switch (Metrics[1],
                         rmse = which.min(Obj),
                         cor = which.max(Obj))
      Result <- genomewideregression(Y,
                                     ETA,
                                     Priortype,
                                     Methodcode,
                                     P,
                                     Hyperparameters[Bestset,],
                                     Thresholdvalue,
                                     Maxiteration,
                                     Randomize,
                                     Variance,
                                     Division_H,
                                     Division_V)
    }
    if(Nset.all > 1){
      Result$Metrics <- cbind(1:Nset.all, Obj, Hyperparameters)
      Colnames <- switch (Metrics[1],
                          rmse = c("Set", "RMSE"),
                          cor = c("Set", "Cor"))
      for(method in 1:Nm)
        Colnames <-
        c(Colnames,
          switch (ETA[[method]]$model,
                  BL     = c("BL_Phi", "BL_Omega"),
                  EBL    = c("EBL_Phi", "EBL_Omega", "EBL_Psi", "EBL_Theta"),
                  BayesA = c("BayesA_Nu", "BayesA_S2", "BayesA_Kappa"),
                  BayesB = c("BayesB_Nu", "BayesB_S2", "BayesB_Kappa"),
                  BayesC = c("BayesC_Nu", "BayesC_S2", "BayesC_Kappa"),
                  BRR    = c("BRR_Nu", "BRR_S2", "BRR_Kappa"),
                  BLUP   = c("BLUP_Nu", "BLUP_S2"))
        )
      colnames(Result$Metrics) <- Colnames
      rownames(Result$Metrics) <- NULL
      Result$Metrics <- data.frame(Result$Metrics)
    }
  }#tuning

  if(Function[1] == "cv"){
    if(Verbose) cat("Cross-validation\n")
    Prediction <- numeric(Nall)
    Obj <- NULL
    for(fold in 1:Nfold){
      if(Verbose){cat("CV fold", fold, "\n")}
      Test <- sort(Partition[, fold])
      Test <- Test[Test != -9]

      ETA.fold <- ETA
      for(method in 1:Nm){
        if(ETA.fold[[method]]$model == "BLUP"){
          ETA.fold[[method]]$K <- ETA[[method]]$K[-Test, -Test, drop = FALSE]
          K.eigen <- eigen(ETA.fold[[method]]$K)
          ETA.fold[[method]]$iK <-
            K.eigen$vectors %*% diag(1/K.eigen$values) %*% t(K.eigen$vectors)
          ETA.fold[[method]]$values <- 1/K.eigen$values
          ETA.fold[[method]]$vectors <- K.eigen$vectors
          ETA.fold[[method]]$tvectors <- t(K.eigen$vectors)
        }else{
          ETA.fold[[method]]$X <- ETA.fold[[method]]$X[-Test, , drop = FALSE]
        }
      }

      if(Nset.all == 1){
        Result <- genomewideregression(Y[-Test],
                                       ETA.fold,
                                       Priortype,
                                       Methodcode,
                                       P,
                                       Hyperparameters,
                                       Thresholdvalue,
                                       Maxiteration,
                                       Randomize,
                                       Variance,
                                       Division_H,
                                       Division_V)
      }else{
        Obj.fold <- tuning(Y[-Test],
                           ETA.fold,
                           Priortype,
                           Methodcode,
                           P,
                           Hyperparameters,
                           Thresholdvalue,
                           Maxiteration,
                           Randomize,
                           CVFoldTuning,
                           Variance,
                           Division_H,
                           Division_V,
                           Metrics[1],
                           Verbose)
        Bestset <- switch (Metrics[1],
                           rmse = which.min(Obj.fold),
                           cor = which.max(Obj.fold))
        Result <- genomewideregression(Y[-Test],
                                       ETA.fold,
                                       Priortype,
                                       Methodcode,
                                       P,
                                       Hyperparameters[Bestset,],
                                       Thresholdvalue,
                                       Maxiteration,
                                       Randomize,
                                       Variance,
                                       Division_H,
                                       Division_V)
      }

      for(method in 1:Nm){
        if(Methodcode[method] == 8){
          #BLUP
          Prediction[Use][Test] <-
            Prediction[Use][Test] +
            ETA[[method]]$K[Test, -Test] %*%
            ETA.fold[[method]]$iK %*%
            matrix(Result$ETA[[method]]$U, ncol = 1)
        }else{
          Prediction[Use][Test] <-
            Prediction[Use][Test] +
            ETA[[method]]$X[Test, , drop=FALSE] %*%
            matrix(Result$ETA[[method]]$Beta, ncol = 1)
        }
      }

      if(Nset.all > 1){
        Obj <- rbind(Obj,
                     c(fold,
                       Bestset,
                       Obj.fold[Bestset],
                       Hyperparameters[Bestset, ]))
      }
    }#fold

    if(Nset.all > 1){
      Colnames <- switch (Metrics[1],
                          rmse = c("Fold", "Set", "RMSE"),
                          cor = c("Fold", "Set", "Cor"))
      for(method in 1:Nm)
        Colnames <-
          c(Colnames,
            switch (ETA[[method]]$model,
                    BL     = c("BL_Phi", "BL_Omega"),
                    EBL    = c("EBL_Phi", "EBL_Omega", "EBL_Psi", "EBL_Theta"),
                    BayesA = c("BayesA_Nu", "BayesA_S2", "BayesA_Kappa"),
                    BayesB = c("BayesB_Nu", "BayesB_S2", "BayesB_Kappa"),
                    BayesC = c("BayesC_Nu", "BayesC_S2", "BayesC_Kappa"),
                    BRR    = c("BRR_Nu", "BRR_S2", "BRR_Kappa"),
                    BLUP   = c("BLUP_Nu", "BLUP_S2"))
        )
      colnames(Obj) <- Colnames
      rownames(Obj) <- NULL
      Obj <- data.frame(Obj)
      if(Random) {
        Result <- list(Prediction = Prediction,
                       Metrics = Obj,
                       Partition = Partition)
      }else{
        Result <- list(Prediction = Prediction, Metrics = Obj)
      }
    } else {
      if(Random) {
        Result <- list(Prediction = Prediction, Partition = Partition)
      }else{
        Result <- list(Prediction = Prediction)
      }
    }
  }#cv

  Result$AddIntercept <- AddIntercept

  if(Verbose) cat("Finished\n")
  Result
}
