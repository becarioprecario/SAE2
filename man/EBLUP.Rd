\name{EBLUP}

\alias{EBLUP}
\alias{EBLUP.ML}
\alias{EBLUP.REML}
\alias{SEBLUP}
\alias{SEBLUP.ML}
\alias{SEBLUP.REML}

\title{EBLUP estimator for small areas}

\description{

This function compute the EBLUP estimator for small areas (Rao, 2003) using a
area level model (see references below). Either Maximum Likelihood (ML) or
Restricted (REML) can be used in the estimation.

You can find more information and examples in the vignette included in the
package.
}



\usage{
EBLUP(formula, varformula, data, Z=NULL, tol=10e-5, maxiter=50, method="ML", na.action=NULL)
SEBLUP(formula, varformula, data, Z=NULL, W, tol=10e-5, maxiter=50, method="ML", na.action=NULL)
}


\arguments{
\item{formula}{Formula relating the response to the covariates in the fixed effects.}
\item{varformula}{A one sided formula indicating the variable which contains thesampling variance.}
\item{data}{A data frame containing the variables used in the model.}
\item{Z}{Structure of the random effects. By default, the identity matrix is used.}
\item{W}{Adjacency matrix used in the Spatial EBLUP.}
\item{tol}{Tolerance used in the computations.}
\item{maxiter}{Maximum number of iterations of the fitting algorithm.}
\item{method}{Either "ML" or "REML"}.
\item{na.action}{Action to handle missing values. At the moment, nothing is 
implemented.}
}

\value{
\item{EBLUP}{EBLUP estimates of the area values.}
\item{beta}{Estimates of the coefficients of the fixed part of the model.}
\item{sigma2u}{Estimate of the variance of the random effects.}
\item{g1}{Estimate of first component (G1) of the MSE.}
\item{g2}{Estimate of second component (G2) of the MSE.}
\item{g3}{Estimate of third component (G3) of the MSE.}
\item{mse}{Estimate of the MSE.}
\item{randeff}{Estimates of the random effects.}
\item{varbeta}{Estimate of the variance of the coefficients of the fixed part.}
}


\seealso{
SEBLUP
}

\references{
JNK Rao (2003).\emph{Small Area Estimation}. John Wiley & Sons, Inc., Hoboken, 
New Jersey.
}

\examples{
#Load data
data(seblup)

spam.options(eps=.0000000001)

d<-data.frame(direst=ydir1, covariate=Xpop[,2], desvar=vardir)

eblupml<-EBLUP(direst~covariate, ~desvar, d)
eblupreml<-EBLUP(direst~covariate, ~desvar, d, method="REML")

summary(eblupreml$randeff - eblupml$randeff)

plot(eblupreml$randeff , eblupml$randeff)
abline(0,1)


plot(eblupml$mse, eblupreml$mse)
abline(0,1)

#Spatial EBLUP
seblupml<-SEBLUP(direst~covariate, ~desvar, d, W=W)
seblupreml<-SEBLUP(direst~covariate, ~desvar, d, W=W, method="REML")

}


\keyword{models}

