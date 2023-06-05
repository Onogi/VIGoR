\name{vigor}
\alias{vigor}
\title{Variational Bayesian inference for genome-wide regression}
\description{
This function performs Bayesian genome-wide regression using variational Bayesian algorithms.
The available regression methods are Bayesian lasso (BL), extended Bayesian lasso (EBL),
BayesA, BayesB, BayesC, Bayesian ridge regression (BRR), BLUP, and fixed effects (FIXED) (fixed effects mean regression using noninformative priors). This function allows multiple regression methods (learners) in a single model. For example, additive effects and interaction effects can be incorporated in a single model using BL and BayesB which provide different shrinkage levels.
}
\usage{
vigor(Y, ETA, Function = c("fitting", "tuning", "cv"), Nfold = 5, CVFoldTuning = 5,
      Partition = NULL, Thresholdvalue = 1e-5,
      Maxiteration = 1000, RandomIni = TRUE, Metrics = c("rmse", "cor"), Verbose = TRUE)
}
\arguments{
  \item{Y}{
  An N-length vector of response variables, where N is the number of samples. Missing data (coded as NA) are allowed.
  }
  \item{ETA}{
  A nested list to specify regression methods, explanatory variables, and hyperparameters. The length of ETA is the number of methods (learners) incorporated in a single model. See details below.
  }
  \item{Function}{
  One of the strings "fitting", "tuning", and "cv". See details below.
  }
  \item{Nfold}{
  An integer value. When n > 1, n-fold cross-validation (CV) is performed on randomly partitioned individuals. When the integer is -1, leave-one-out CV is conducted. Used when Function = "cv" and Partition == NULL.
  }
  \item{CVFoldTuning}{
  An integer specifying the fold number of the CV in hyperparameter tuning. Used when Function = "cv" or "tuning" and multiple hyperparameter sets are given.
  }
  \item{Partition}{
  A matrix defining the partitions of CV. See details and examples below. Used when Function = "cv".
  }
  \item{Thresholdvalue}{
  Specifies the convergence threshold. Smaller values indicate stricter thresholds.
  }
  \item{Maxiteration}{
  Maximum number of iterations.
  }
  \item{RandomIni}{
  If TRUE, the initial values of the SNP effects are randomly determined. Otherwise, they are set to 0.
  }
  \item{Metrics}{
  One of the strings "rmse" and "cor" to specify the metrics used in CV. rmse and cor use RMSE and Pearson correlation, respectively.
  }
  \item{Verbose}{
  If TRUE, print the run information to the console.
  }
}
\details{
  \bold{Regression methods}\cr
  Vigor supports the following regression methods;
  \itemize{
    \item BL (Bayesian lasso)
    \item EBL (extended Bayesian lasso)
    \item BayesA
    \item BayesB
    \item BayesC
    \item BRR (Bayesian ridge regression)
    \item BLUP
    \item FIXED (fixed effects)
  }
  These methods can be included in a single model simultaneously with different explanatory variables. For the details of these methods and the theoretical backgrounds of vigor, see the pdf document (Onogi 2021). \cr

  \bold{ETA}\cr
  Each element (list) of ETA consists of the following objects.\cr
    \itemize{
    \item model : One of strings representing regression methods, "BL", "EBL", "BayesA", "BayesB", "BayesC", "BRR", "BLUP", or "FIXED"
    \item X     : An explanatory variables (e.g., SNP genotypes) of (N x P) matrix, where N and P denote the number of samples and variables, respectively
    \item K     : An N x N kernel matrix (e.g., genomic relationship matrix) used when BLUP is specified as model. When BLUP is specified but K is lacked, the linear kernel is created from X as scale(X)*t(scale(X))/ncol(X).
    \item H     : A vector or matrix including hyperparameters
    \item data  : A data frame containing fixed effects. Used to model FIXED using formula.
    }
  Specification of model is essential for all methods. For regression methods except for BLUP and FIXED, X is essential. For BLUP, either X or K is essential. For FIXED, either X or formula with data is essential.

  \bold{Hyperparameters}\cr
  The regression methods require hyperparameters as H in ETA. H can be a vector or matrix. Below is the order of hyperparameters in H (order of columns in the case of matrix). Default values are shown in parenthesis.\cr
  \itemize{
    \item BL     : Phi(1), Omega(1)
    \item EBL    : Phi(0.1), Omega(0.1), Psi(1), Theta(0.1)
    \item BayesA : Nu(5), S2(0.01)
    \item BayesB : Nu(5), S2(0.1), Kappa(0.01)
    \item BayesC : Nu(5), S2(0.1), Kappa(0.01)
    \item BRR    : Nu(5), S2(0.01)
    \item BLUP   : Nu(5), S2(0.3)
    }
  \bold{Note that Kappa is the proportion of explanatory variables with NON-ZERO EFFECTS}. Also note that Y is standardized automatically. To specify multiple hyperparameter sets, give an S x Nh matrix where S is the number of sets and Nh is the number of hyperparameters of the method to ETA. See the pdf document (Onogi 2021) for the details of hyperparameters.\cr

 \bold{Functions}\cr
  The functions of vigor are "fitting", "tuning", and "cv".
    \itemize{
      \item "fitting" : Fits the specified regression model to the data. When H (hyperparameters) includes multiple hyperparameter sets (i.e., H is a matrix), only the first set is used.
      \item "tuning"  : Selects the best hyperparameter set using CV for tuning. This set is then used for model fitting. The CV is performed on randomly partitioned data. The number of folds is determined by CVFoldTuning. When multiple hyperparameter sets are given to multiple methods, all combinations of hyperparameters are attempted.
      \item "cv"      : Conducts CV and returns the predicted values. When multiple hyperparameter sets are given, tuning is performed at each fold of the CV. As in Tuning, when multiple hyperparameter sets are given to multiple methods, all combinations of hyperparameters are attempted.
    }

\bold{Partition matrix}\cr
  The following is a possible Partition of 20 individuals evaluated in a five-fold CV:\cr\cr
    14 11  3  2  7\cr
     5  4 20 10  9\cr
     6  8 16 15 12\cr
    18 13 17  1 19\cr\cr
  Sample (row numbers in Y/X/K) 14, 5, 6, and 18 are removed from the training set at the first fold of the five-fold CV. Samples 11, 4, 8, and 13 are removed at the next fold. This process is repeated up to the fold number of the CV. If the number of samples N is 19, the gap is filled with -9. For example,\cr\cr
     8  6  3 14 18\cr
    12  4  1 15  5\cr
    17  9 13 11 10\cr
    19 16  7  2 -9\cr\cr
  An example of random sampling validation in which samples can be sampled more than once is shown below.\cr\cr
    18  3 11 16 13\cr
    17  8 13 13 18\cr
     7 15 14 19  7\cr
     1 13 12  7  2\cr\cr
  Samples 18, 13, and 7 are repeatedly used as testing samples.\cr\cr
  Random partitioning outputs a Partition matrix, which can be input as the Partition matrix in subsequent analysis.\cr

  \bold{Intercept}\cr
  If No FIXED is given in ETA, vigor automatically adds the intercept to the regression model as a fixed effect (FIXED). If FIXED is given by the user, vigor regards the first column as the intercept\cr

  \bold{Standardization}\cr
  Vigor standardizes Y (response variables). Although most returned values are scaled back to the original scale, the lower bound of the marginal log likelihood is returned as the standardized scale.
}

\value{
  When Function = "fitting" or "tuning", a list containing the following elements is returned.
  \item{$LB}{Lower bound of the marginal log likelihood of Y.}
  \item{$ResidualVar}{Residual variances (1/Tau02) at each iteration (from start to end).}
  \item{$H}{Used hyperparameters.}
  \item{$Fittedvalue}{Fitted values.}
  \item{$Metrics}{Metrics for hyperparameter tuning. Returned when Function = "tuning"}
  \item{$ETA}{A list containing results for each method.
    \itemize{
    \item Beta: Posterior means of regression coefficients of X
    \item Sd.beta: Posterior standard deviations of Beta (uncertainty of Beta)
    \item Sigma2: Posterior means of variance of Beta or U
    \item Rho: Posterior means of model-inclusion probabilities
    \item U: Posterior means of random effects of BLUP
    \item Sd.u: Posterior standard deviations of U (uncertainty of U)
    \item iK: Inverse of K
    }
    Beta and Sd.beta are returned for BL, EBL, BayesA, BayesB, BayesC, and FIXED, and U, Sd.u, and iK are returned for BLUP. Rho is returned for BayesB and BayesC. Sigma2 is returned for methods except for BL and EBL.}
  \item{$AddIntercept}{True when the intercept was added automatically.}
  When Function = "cv", a list containing the following elements is returned.
  \item{$Prediction}{A vector of predicted values}
  \item{$Metrics}{Metrics of hyperparameter tuning. Chosen sets and corresponding metrics at each fold are returned.}
  \item{$Partition}{A matrix representing the partition used in random partitioning. This matrix can be used as the argument Partition in subsequent analyses.}
  \item{$AddIntercept}{True when the intercept was added automatically.}
}
\references{
  Onogi A., Variational Bayesian inference for genome-wide regression: joint estimation of multiple learners, in prep.\cr
  Onogi A. & Iwata H., 2016 VIGoR: Variational Bayesian Inference for Genome-Wide Regression.
  Journal of Open Research Software, 4: e11\cr
  Onogi A., 2021, Documents for VIGoR ver. 1.1.0, https://github.com/Onogi/VIGoR\cr
}
\author{
  Akio Onogi
}
\examples{
\donttest{
#DATA###########################################################################
data(sampledata)
dim(X) #Matrix of SNP genotypes (explanatory variables)
dim(Z) #Matrix of a fixed effect (explanatory variables)
length(Y) #Vector of response variables


#Fitting########################################################################
#Example 1: Fit SNP genotypes with BayesC
ETA <- list(list(model = "BayesC", X = X))
Result <- vigor(Y, ETA)
##see estimated SNP effects
plot(abs(Result$ETA[[1]]$Beta), pch = 20)
##10 SNPs at 1, 101, ..., 901 have non-zero effects
abline(v = seq(1, 1000, 100), col = 2, lty = 2)
##see inclusion probability
plot(Result$ETA[[1]]$Rho, pch = 20)
abline(v = seq(1, 1000, 100), col = 2, lty = 2)
##Intercept is added automatically as the last learner
Result$ETA[[2]]$Beta


#Example 2: Fit fixed effects and SNP genotypes
##There are two approaches to fit fixed effects
##(1) Create model matrix
Z #Z consists of three categories (A, B, and C)
Z.matrix <- model.matrix(~ Z)
head(Z.matrix) #The first column is the intercept
##Fit with EBL
ETA <- list(list(model = "FIXED", X = Z.matrix),
            list(model = "EBL", X = X))
Result <- vigor(Y, ETA)
##Estimated fixed effects (intercept, B, and C)
Result$ETA[[1]]$Beta
##Estimated SNP effects
plot(abs(Result$ETA[[2]]$Beta), pch = 20)
abline(v = seq(1, 1000, 100), col = 2, lty = 2)
##NOTE: when FIXED is added by user, the intercept is not automatically added.
##Thus, variables in FIXED should contain the intercept.

##(2) Use formula
Data <- data.frame(Z = factor(Z))
ETA <- list(list(~ Z, model = "FIXED", data = Data),
            list(model = "BayesA", X = X))
Result <- vigor(Y, ETA)
##Estimated fixed effects (intercept, B, and C)
Result$ETA[[1]]$Beta
plot(abs(Result$ETA[[2]]$Beta), pch = 20)
abline(v=seq(1,1000,100),col=2,lty=2)
##NOTE: formula automatically adds the intercept


#Example 3: Multiple regression methods in a single model
##Some SNPs in X have dominance (non-additive) effects
##Fit SNP genotypes coded as additive and dominance with different shrinkage levels
X.d <- X
X.d[X == 2] <- 0 #heterozygotes are 1 and homozygotes are 0
ETA <- list(list(~ Z, model = "FIXED", data = Data),
            list(model = "BayesC", X = X, H = c(5, 0.1, 0.01)),
            list(model = "BayesC", X = X.d, H = c(5, 0.1, 0.001)))
Result <- vigor(Y, ETA)
##Inclusion probability for additive effects
plot(Result$ETA[[2]]$Rho, pch = 20)
abline(v = seq(1, 1000, 100), col = 2, lty = 2)
##Inclusion probability for dominance effects
plot(Result$ETA[[3]]$Rho, pch = 20)
##SNPs at 1, 201, ..., 801 have non-zero effects
abline(v = seq(1, 1000, 200), col = 2, lty = 2)

##Fit additive and dominance effects with different learners
ETA <- list(list(~ Z, model = "FIXED", data = Data),
            list(model = "BL", X = X, H = c(1, 0.01)),
            list(model = "BayesC", X = X.d, H = c(5, 0.1, 0.001)))
Result <- vigor(Y, ETA)
plot(abs(Result$ETA[[2]]$Beta), pch = 20)
abline(v = seq(1, 1000, 100), col = 2, lty = 2)
plot(Result$ETA[[3]]$Rho, pch = 20)
abline(v = seq(1, 1000, 200), col = 2, lty = 2)


#Tuning hyperparameters#########################################################
#Example 4: Model fitting after hyperparameter tuning with cross-validation
##Candidate hyperparameter values are determined with hyperpara
##Use BayesB
ETA <- list(list(~ Z, model = "FIXED", data = Data),
            list(model = "BayesB", X = X,
                 H = hyperpara(X, 0.5, "BayesB", c(0.1,0.01))))
Result <- vigor(Y, ETA, Function = "tuning")
##See the tuned result
Result$Metrics
##The model was fitted to the full data with the best set
Result$H
plot(Result$ETA[[2]]$Rho, pch = 20)
abline(v = seq(1, 1000, 100), col = 2, lty = 2)
##See how much the model was fitted
plot(Y, Result$Fittedvalue); abline(0, 1)

##When multiple learners used, all combinations of hyperparameter sets are compared
ETA <- list(list(~ Z, model = "FIXED", data = Data),
            list(model = "BayesB", X = X,
                 H = hyperpara(X, 0.5, "BayesB", c(0.1,0.01))),
            list(model = "BayesC", X = X.d,
                 H = hyperpara(X, 0.5, "BayesC", c(0.1,0.01))))
Result <- vigor(Y, ETA, Function = "tuning")
##BayesB and BayesC have two candidate sets, respectively.
##Thus, total 2 x 2 = 4 combinations are compared.
Result$Metrics
##The model was fitted to the full data with the best combination.
Result$H
plot(Result$ETA[[2]]$Rho, pch = 20)
abline(v = seq(1, 1000, 100), col = 2, lty = 2)
plot(Result$ETA[[3]]$Rho, pch = 20)
abline(v = seq(1, 1000, 200), col = 2, lty = 2)


#Cross-validation###############################################################
#Example 5: Cross-validation with random splitting
ETA <- list(list(~ Z, model = "FIXED", data = Data),
            list(model = "BayesC", X = X,
                 H = hyperpara(X, 0.5, "BayesC", c(0.1,0.01))))
Result <- vigor(Y, ETA, Function="cv")
##Because two hyperparameter sets are provided,
##nested CV is conducted at each fold to tune hyperparameters
##See which set was selected at each fold
Result$Metrics
##See predicted values
plot(Y, Result$Prediction)
cor(Y, Result$Prediction)

#Example 6: Cross-validation with the specified splitting
##Perform CV using the same partition as Example 5. Use EBL
ETA <- list(list(~ Z, model = "FIXED", data = Data),
            list(model = "EBL", X = X))
Result2 <- vigor(Y, ETA, Function="cv", Partition = Result$Partition)
plot(Y, Result2$Prediction)
cor(Y, Result2$Prediction)
}
}
