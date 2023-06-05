\name{predict_vigor}
\alias{predict_vigor}
\title{
Predict Y of new data using a training result with vigor
}
\description{
This function predicts Y of new data using a training result with vigor.
}
\usage{
predict_vigor(training_result, newX)
}
\arguments{
  \item{training_result}{
  Result of vigor
}
  \item{newX}{
  A list of X. newX contains X for each learner. The length and order of learners should be same as those of ETA used for training. Each element of newX is a matrix. When BLUP is used, X is an n1 x n2 relationship matrix where n1 and n2 are the numbers of samples in test and training data, respectively. When the other methods are used, X is an n1 x p matrix where p is the number of explanatory variables. p should be the same as the training data. When the intercept is added automatically, newX needs not to include the intercept.
}
}
\details{
This function predict Y of new data (newX). When multiple learners are included in the model, predicted values are calculated for each element of newX, and the summation of all predicted values is returned (predicted values of each learner are not returned). When contributions of each learner (method) are of interest, add NULL to the elements of newX to be ignored (see examples).
}
\value{
A vector of predicted values is returned.
}
\seealso{
vigor
}
\examples{
\donttest{
#data
data(sampledata)
dim(X) #Matrix of SNP genotypes (explanatory variables)
dim(Z) #Matrix of a fixed effect (explanatory variables)
length(Y) #Vector of response variables

#Train EBL using the first 80 percent of data
#Then predict the remaining 20 percent data
ETA <- list(list(model = "EBL",X = X[1:400, ]))
Train <- vigor(Y[1:400], ETA)
newX <- list(X[401:500, ])
Predict <- predict_vigor(Train, newX)
plot(Y[401:500], Predict)
#When the intercept is automatically added when training,
#the intercept is again automatically added to predicted values

#Use multiple regression methods
#Fit additive and dominance effects using BayesC with different shrinkage levels
#Also fixed effects are added
X.d <- X
X.d[X == 2] <- 0 #heterozygotes are 1 and homozygotes are 0
Z.matrix <- model.matrix(~ Z)
ETA <- list(list(model = "FIXED", X = Z.matrix[1:400, ]),
            list(model = "BayesC", X = X[1:400, ], H = c(5, 0.1, 0.01)),
            list(model = "BayesC", X = X.d[1:400, ], H = c(5, 0.1, 0.001)))
Train <- vigor(Y[1:400], ETA)
newX <- list(Z.matrix[401:500, ], X[401:500, ], X.d[401:500, ])
Predict <- predict_vigor(Train, newX)
plot(Y[401:500], Predict)

#When fixed effects are specified using formula,
Data <- data.frame(Z = factor(Z))
ETA <- list(list(~ Z, model = "FIXED", data = Data[1:400, ,drop=FALSE]),
            list(model = "BayesB", X = X[1:400, ], H = c(5, 0.1, 0.01)),
            list(model = "BayesB", X = X.d[1:400, ], H = c(5, 0.1, 0.001)))
Train <- vigor(Y[1:400], ETA)
newX <- list(~ Z, data = Data[401:500, , drop=FALSE], X[401:500, ], X.d[401:500, ])
Predict <- predict_vigor(Train, newX)
plot(Y[401:500], Predict)
#NOTE: please confirm that levels of fixed effects are consistent
#between the training and testing data

#Contributions of each learner can be assessed by filling newX with NULL
##Contribution of additive effect
newX <- list(NULL, X[401:500, ], NULL)
Predict <- predict_vigor(Train, newX)
plot(Y[401:500], Predict)

##Contribution of dominance effect
newX <- list(NULL, NULL, X.d[401:500, ])
Predict <- predict_vigor(Train, newX)
plot(Y[401:500], Predict)
}
}
