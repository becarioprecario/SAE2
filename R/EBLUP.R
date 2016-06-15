
EBLUP <- function(formula, varformula, data, Z=NULL, tol=10e-5, maxiter=50, method="ML", na.action=NULL)
{

	mf<-model.frame(formula, data)

	direst<-matrix(model.response(mf), ncol=1)
	X<-model.matrix(formula, data)
	desvar<-matrix(model.matrix(varformula, data)[,2], ncol=1)
	k<-nrow(direst)

	if(is.null(Z))
		Z<-diag.spam(1,k)


	res<-switch(method,
		ML = EBLUP.ML(direst, X, desvar, tol=tol, maxiter=maxiter),
		REML = EBLUP.REML(direst, X, desvar, tol=tol, maxiter=maxiter)
	)

	if(is.null(res))
		print("Method should be ML or REML\n")

	#Add some other stuff to the object

	res$coef<-as.matrix(res$coef)
	dimnames(res$coef)<-list(dimnames(X)[[2]], NULL)

	res$rank<-min(dim(X))
	#res$df.residual<-nrow(direst)-res$rank#Is this right?
	res$call<-match.call()

	res$Z<-Z
	res$desvar<-desvar

	#We may want to provide a better way of handling missing data
	res$na.action<-na.action	

	return(res)


}

EBLUP.ML<-function(direst, X, desvar, Z=NULL, tol=10e-5, maxiter=50)
{

	k<-nrow(direst)

	R<-diag.spam(desvar[,1], k)

	if(is.null(Z))
		Z=diag.spam(1, k)

	l<- t(X) #FIXME!! Make sure that this l what it should be
	m<- Z    #FIXME!! Make sure that this l what it should be

	#Iterative algorithm to estimate variances goes here

	niter<-0
	dif<-tol+1

	ZtZ<-Z%*%t(Z)


	#Variance of the residuals as initial point	
	sigma2u<-as.numeric(var(direst -X%*%solve(t(X)%*%diag(1/desvar[,1])%*%X)%*%t(X)%*%diag(1/desvar[,1])%*%direst))

	while((niter<maxiter) & (dif > tol))
	{
		G<-diag.spam(sigma2u, k)
		GtZ<-G%*%t(Z)
		V<-R + Z%*%GtZ

		Vinv<-as.spam(solve(V))
		VinvX<-Vinv%*%X
		slvXVinvX<-solve(t(X)%*%VinvX)

		coeff<-slvXVinvX%*%t(VinvX)%*%direst

		VinvZtZ<-Vinv%*%ZtZ
		#Partial derivative of log-likelihood of sigma2u
		s<- (-.5)*sum(diag(VinvZtZ)) -.5*t(direst-X%*%coeff)%*%(-VinvZtZ%*%Vinv) %*%(direst-X%*%coeff)

		#Matrix of expected second derivaties of sigma2u
		Isigma2u<-.5*sum(diag(VinvZtZ%*%VinvZtZ))

		sigma2u.old<-sigma2u
		sigma2u<-as.numeric(sigma2u+ (1/Isigma2u)*s)

		niter<-niter+1
		dif<-abs(sigma2u.old-sigma2u)

	}
	sigma2u<-sigma2u.old

#	G<-diag.spam(sigma2u, k)
#	GtZ<-G%*%t(Z)
#	V<-R + Z%*%GtZ
#
#	Vinv<-as.spam(solve(V))
#	VinvX<-Vinv%*%X
#
#	slvXVinvX<-solve(t(X)%*%VinvX)
#
#	coeff<-slvXVinvX%*%t(VinvX)%*%direst

	ftd<- X%*%coeff

	randeff<-G%*%t(Z)%*%Vinv%*%(direst-ftd)


#	Compute g1, g2 and g3 here

	g1<- diag(t(m)%*%(G-GtZ%*%Vinv%*%t(GtZ))%*%m)

	#Bias correction for g1. As in Rao, 2003, page 109
	Vinvder<- (-VinvZtZ%*%Vinv)
	bsigma2u<-(1/(2*k))*(1/Isigma2u)*sum(diag( slvXVinvX %*% (t(X)%*%Vinvder%*%X ) ))
	g1grad<-diag(  t(m)%*%m - t(m)%*%( t(Z)%*%Vinv%*%t(GtZ) - GtZ%*%Vinvder%*%t(GtZ) +GtZ%*%Vinv%*%Z )%*%m  )
	g1biascorr<- bsigma2u*g1grad


	g2<-sapply(1:k, function(i)
		{
			tb<-t(m[,i])%*%GtZ%*%VinvX
			td<-t(l[,i])-tb

			td%*%slvXVinvX%*%t(td)
		})

	VinvZ<-Vinv%*%Z

	g3<-sapply(1:k, function(i)
		{

	deltab<-t(m[,i])%*%( t(VinvZ) + GtZ%*%(-VinvZ%*%t(VinvZ)) )

	sum(diag(deltab%*%V%*%t(deltab))*(1/Isigma2u)) 

		})


	res<-list(
		coefficients=coeff, 
		residuals=as.vector((direst-ftd)[,1]),
		#effects =NULL,
		#rank=NULL,
		fitted.values=as.vector(ftd[,1]),
		randeff=as.vector(randeff[,1]),
		varcoeff=slvXVinvX,
		varsigma2u=1/Isigma2u,
		sigma2u=sigma2u,
		g1=g1,
		g2=g2,
		g3=g3,
		g1biascorrd=g1biascorr,
		mse= g1-g1biascorr+g2+2*g3,
		eblup=as.vector(ftd[,1]+randeff[,1]),
		method="ML"
	)

	class(res)<-"EBLUP"

	return(res)
}








EBLUP.REML<-function(direst, X, desvar, Z=NULL, tol=10e-5, maxiter=50)
{

	#Variance of the residuals as initial point	
	sigma2u<-as.numeric(var(direst -X%*%solve(t(X)%*%diag(1/desvar[,1])%*%X)%*%t(X)%*%diag(1/desvar[,1])%*%direst))

	k<-nrow(direst)

	R<-diag.spam(desvar[,1], k)

	if(is.null(Z))
		Z=diag.spam(1, k)

	l<- t(X) #FIXME!! Make sure that this l what it should be
	m<- Z    #FIXME!! Make sure that this l what it should be

	#Iterative algorithm to estimate variances goes here

	niter<-0
	dif<-tol+1

	ZtZ<-Z%*%t(Z)
	while((niter<maxiter) & (dif > tol))
	{
		G<-diag.spam(sigma2u, k)
		GtZ<-G%*%t(Z)
		V<-R + Z%*%GtZ

		#V^{-1} and 
		Vinv<-as.spam(solve(V))
		VinvX<-Vinv%*%X
		slvXVinvX<-solve(t(X)%*%VinvX)


		P<-Vinv-VinvX%*%slvXVinvX%*%t(VinvX)

		#coeff<-slvXVinvX%*%t(VinvX)%*%direst

		VinvZtZ<-Vinv%*%ZtZ
		#Partial derivative of log-likelihood of sigma2u
		PZtZ<-P%*%ZtZ
		s<- (-.5)*sum(diag(PZtZ)) +.5*t(direst)%*%(PZtZ%*%P)%*%direst

		#Matrix of expected second derivaties of sigma2u
		Isigma2u<-.5*sum(diag(PZtZ%*%PZtZ))

		sigma2u.old<-sigma2u
		sigma2u<-as.numeric(sigma2u+ (1/Isigma2u)*s)

		niter<-niter+1
		dif<-abs(sigma2u.old-sigma2u)

	}
	sigma2u<-sigma2u.old

#	G<-diag.spam(sigma2u, k)
#	GtZ<-G%*%t(Z)
#	V<-R + Z%*%GtZ
#
#	Vinv<-as.spam(solve(V))
#	VinvX<-Vinv%*%X
#
#	slvXVinvX<-solve(t(X)%*%VinvX)
#
	coeff<-slvXVinvX%*%t(VinvX)%*%direst

	ftd<- X%*%coeff

	randeff<-G%*%t(Z)%*%Vinv%*%(direst-ftd)


#	Compute g1, g2 and g3 here

	g1<- diag(t(m)%*%(G-GtZ%*%Vinv%*%t(GtZ))%*%m)

	#Bias correction for g1. As in Rao, 2003, page 109
#	Vinvder<- (-VinvZtZ%*%Vinv)
#	bsigma2u<-(1/(2*k))*(1/Isigma2u)*sum(diag( slvXVinvX %*% (t(X)%*%Vinvder%*%X ) ))
#	g1grad<-diag(  t(m)%*%m - t(m)%*%( t(Z)%*%Vinv%*%t(GtZ) - GtZ%*%Vinvder%*%t(GtZ) +GtZ%*%Vinv%*%Z )%*%m  )
#	g1biascorr<- bsigma2u*g1grad
	g1biascorr<-0 #FIXME

	g2<-sapply(1:k, function(i)
		{
			tb<-t(m[,i])%*%GtZ%*%VinvX
			td<-t(l[,i])-tb

			td%*%slvXVinvX%*%t(td)
		})

	VinvZ<-Vinv%*%Z

	g3<-sapply(1:k, function(i)
		{

	deltab<-t(m[,i])%*%( t(VinvZ) + GtZ%*%(-VinvZ%*%t(VinvZ)) )

	sum(diag(deltab%*%V%*%t(deltab))*(1/Isigma2u)) 

		})


	res<-list(
		coefficients=coeff, 
		residuals=as.vector((direst-ftd)[,1]),
		#effects =NULL,
		#rank=NULL,
		fitted.values=as.vector(ftd[,1]),
		randeff=as.vector(randeff[,1]),
		varcoeff=slvXVinvX,
		sigma2u=sigma2u,
		varsigma2u=1/Isigma2u,
		g1=g1,
		g2=g2,
		g3=g3,
		g1biascorrd=g1biascorr,
		mse= g1-g1biascorr+g2+2*g3,
		eblup=as.vector(ftd[,1]+randeff[,1]),
		method="REML"
	)
	class(res)<-"EBLUP"

	return(res)
}




residuals.EBLUP <- function(object, ...) {
        if (is.null(object$na.action))
                object$residuals
        else napredict(object$na.action, object$residuals)
}

fitted.EBLUP <- function(object, ...) {
        if (is.null(object$na.action))
                object$fitted.values
        else napredict(object$na.action, object$fitted.values)
}
coef.EBLUP <- function(object, ...) {
        object$coefficients
}

print.EBLUP <- function(x, ...) {
        cat("\nCall:\n")
        print(x$call)
        cat("\nCoefficients:\n")
        print(coef(x))
	cat("\nVariance of the random effects:", x$sigma2u, "\n")
        cat("\nLog likelihood:", logLik(x), "\n")
        invisible(x)

}


ranef.EBLUP <-function(object, ...) {
	object$randeff
}


#FIXME: What likelihood shall we return? Profile, conditional, etc.?
logLik.EBLUP <- function(object, ..., conditional=TRUE) {


	if(conditional)
	{
	   mres<-residuals(object) - object$Z%*%ranef(object)
	   V<-diag.spam(object$desvar[,1])
	}		
	else
	{
		mres<-matrix(fitted(object), ncol=1)
		V<-diag.spam(object$desvar[,1])+object$Z%*%diag(object$sigma2u, nrow=nrow(object$desvar))%*%t(object$Z)
	}
	
	n<-nrow(mres)
	l<- (-n/2)*determinant(V)$modulus -.5*t(mres)%*%solve(V)%*%mres	
	
	return(l[1,1])
}

deviance.EBLUP <- function(object, ...) {
       return(-2*logLik(object)) 
}

AIC.EBLUP<-function(object, ..., k=2) {
	return(deviance(object, ...)+k*(object$rank+1) )#+1 because of var. r. effects
}


#FIXME: Implement cAIC
cAIC.EBLUP<-function(object, ..., k=2) {
	return(NA)
}


