\name{hyperpara}
\alias{hyperpara}
\title{
Calculation of hyperparameter values
}
\description{
This function determines the hyperparameter values of regression methods,
based on two assumptions.
}
\usage{
hyperpara(X, Mvar, Model = c("BL", "EBL", "BayesA", "BayesB", "BayesC", "BRR"),
          Kappa = 0.01, Xtype = c("Geno", "Var"), f = 0, BL.Phi = 1, EBL.Phi = 0.1,
          EBL.Omega = 0.1, EBL.Psi = 1, Nu = 5, Verbose=FALSE)
}
\arguments{
  \item{X}{
  An (N x P) matrix, where N and P denote the number of samples and variables, respectively. When X is SNP genotypes as specified by Xtype = "Geno", SNP genotypes should be coded with values between 0 and 2. SNP genotypes are used to calculate sum (2*Q*(1-Q)*(1+f)) where Q is a vector of allele frequencies and f is the inbreeding coefficient. Missing values in X are not allowed. When X is "Var"(variables other than SNPs), variances of variables are used instead of sum (2*Q*(1-Q)*(1+f)).
}
  \item{Mvar}{
  A scalar or vector denoting the assumed proportion of variance of Y that can be explained by X. Mvar is < 1.0 in BL and EBL, and <= 1.0 in the other methods.
}
  \item{Model}{
  One of the six regression methods (BL, EBL, BayesA, BayesB, BayesC, BRR).
}
  \item{Kappa}{
  A scalar or vector containing the assumed proportion of variables with NON-ZERO EFFECTS. Used when BL, EBL, BayesB, and BayesC. Kappa is set to 1 when BayesA and BRR are used.
}
  \item{Xtype}{
  Allowed Xtypes are "Geno" and "Var". Enter "Geno" when X contains the SNP genotypes and "Var" when X contains variables other than SNP genotypes.
}
  \item{f}{
  A scalar representing the inbreeding coefficient. Enter 1 for inbred species.
}
  \item{BL.Phi}{
  A scalar or vector containing Phi values of BL.
}
  \item{EBL.Phi}{
  A scalar or vector containing Phi values of EBL.
}
  \item{EBL.Omega}{
  A scalar or vector containing Omega values of EBL.
}
  \item{EBL.Psi}{
  A scalar or vector containing Psi values of EBL.
}
  \item{Nu}{
  A scalar or vector containing Nu values of BayesA, BayesB, BayesC, and BRR.
}
  \item{Verbose}{
  Specifies whether to print information (TRUE) or not (FALSE).
}
}
\details{
  To run vigor, users must specify the following hyperparameter values.\cr
  \itemize{
    \item BL: Phi, Omega
    \item EBL: Phi, Omega, Psi, Theta
    \item BayesA: Nu, S2
    \item BayesB: Nu, S2, Kappa
    \item BayesC: Nu, S2, Kappa
    \item BRR: Nu, S2
    }
  This function calculates the Omega of BL; Theta of EBL; S2 of BayesA, BayesB, BayesC, and BRR. Mvar, Kappa, and the other hyperparameters required by each method are specified by the user. The definitions of Mvar and Kappa are intuitively understandable and relatively easy to specify (see Arguments). For the other hyperparameters, the default values are recommended. When the arguments of hyperpara are vectors, all value-combinations are returned as a matrix. The hyperparameters are explained in the details of vigor and in the pdf manual of VIGoR (Onogi 2021).
}
\value{
This function returns a vector when yielding a single hyperparameter set, and a matrix when yielding multiple hyperparameter sets. The rows and columns of the matrix correspond to the sets (value combinations) and the hyperparameters, respectively. See examples below.\cr
}
\seealso{
  vigor
}
\references{
  Onogi A., Variational Bayesian inference for genome-wide regression: joint estimation of multiple learners, in prep.\cr
  Onogi A. & Iwata H., 2016 VIGoR: Variational Bayesian Inference for Genome-Wide Regression.
  Journal of Open Research Software, 4: e11\cr
  Onogi A., 2021, Documents for VIGoR ver. 1.1.0, https://github.com/Onogi/VIGoR\cr
}
\examples{
\donttest{
#data
data(sampledata)
dim(X) #500 samples and 1000 variables
unique(X[1:(100*1000)]) #coded as 0, 1, 2

#A single Mvar (0.5) and Kappa (0.01) value is assumed for BL.
#A vector is returned.
hyperpara(X, 0.5, "BL", 0.01, Verbose = TRUE)

#Phi is set to 1 as default. To change Phi, use BL.Phi.
hyperpara(X, 0.5, "BL", 0.01, BL.Phi = 5)

#Calculate multiple hyperparameter value sets of BayesC assuming that Kappa is 0.1 and 0.01.
#A matrix is returned.
hyperpara(X, 0.5, "BayesC", c(0.1, 0.01))

#The output vector can be used as the argument of vigor
ETA <- list(list(model = "BayesB",X = X,
                 H = hyperpara(X, 0.5, "BayesC", c(0.1, 0.01))))
Result <- vigor(Y, ETA, Function = "tuning")
Result$Metrics

#Calculate hyperparameter values of EBL
hyperpara(X, c(0.2, 0.5), "EBL", c(0.1, 0.01), EBL.Omega = c(0.5, 1))
#Total 2 (Mvar) x 2 (Kappa) x 2 (EBL.Omega) = 8 sets are returned.
}
}


