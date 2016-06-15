
SEBLUP <- function(formula, varformula, data, Z=NULL, W, tol=10e-5, maxiter=50, method="ML", na.action=NULL)
{

        mf<-model.frame(formula, data)

        direst<-matrix(model.response(mf), ncol=1)
        X<-model.matrix(formula, data)
        desvar<-matrix(model.matrix(varformula, data)[,2], ncol=1)
        k<-nrow(direst)

        if(is.null(Z))
                Z<-diag.spam(1,k)


        res<-switch(method,
                ML = SEBLUP.ML(direst, X, desvar, Z, W, tol=tol, maxiter=maxiter),
                REML = SEBLUP.REML(direst, X, desvar, Z, W, tol=tol, maxiter=maxiter)
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

	class(res)<-c("SEBLUP", "EBLUP")
        return(res)

}


loglik.SEBLUP<-function(par, direst, X, desvar, Z, W)
{
	ndim<-length(par)
	sigma2u<-par[1]
	rho<-par[2]

	#print(par)	

	k<-nrow(direst)

	V<-diag.spam(desvar[,1])+ sigma2u*Z%*%solve((diag.spam(1,k)-rho*W)%*%(diag.spam(1,k)-rho*t(W)))%*%t(Z)
	Vinv<-as.spam(solve(as.matrix(V)))
	tXVinv<-t(X)%*%Vinv
	coeff<-solve(tXVinv%*%X)%*%tXVinv%*%direst
	l<- ( (-.5)*(determinant(as.matrix(V), logarithm=TRUE))$modulus 
		-.5*t(direst-X%*%coeff)%*%(Vinv)%*%(direst-X%*%coeff) )

	return(-l)

}


gr.SEBLUP<-function(par, direst, X, desvar, Z, W)
{
	sigma2u<-par[1]
	rho<-par[2]


	k<-nrow(X)
	R<-diag.spam(desvar[,1], k)


	s<-c(NA, NA)


	ZtZ<-Z%*%t(Z)
	WtW<-W%*%t(W)

	IrhoW<- diag.spam(1, k)-rho*W
	C<-IrhoW*t(IrhoW)
	slvC<-solve(C)
	G<-sigma2u*C
	GtZ<-G%*%t(Z)
	V<-R + Z%*%GtZ

	Vinv<-as.spam(solve(V))
	VinvX<-Vinv%*%X
	slvXVinvX<-solve(t(X)%*%VinvX)

	coeff<-slvXVinvX%*%t(VinvX)%*%direst

	#VinvZtZ<-Vinv%*%ZtZ

	VZCZ<-Vinv%*%Z%*%slvC%*%t(Z)

	ftd<- X%*%coeff
	rsd<-direst-ftd

	s[1]<-( (-.5)*sum(diag(VZCZ))
		+.5*t(rsd)%*%(VZCZ%*%Vinv)%*%rsd )

	A<-(-sigma2u)*slvC%*%(2*rho*WtW-2*W)%*%(slvC)
	auxs2<-Vinv%*%Z%*%A%*%t(Z)
	s[2]<- ( (-.5)*sum(diag(auxs2))
		+.5*t(rsd)%*%(auxs2%*%Vinv)%*%rsd )

	return(s)
}



SEBLUP.ML<-function(direst, X, desvar, Z=NULL, W, tol=10e-5, maxiter=50)
{

	k<-nrow(direst)

	R<-diag.spam(desvar[,1], k)

	if(is.null(Z))
		Z=diag.spam(1, k)

	l<- t(X)
	m<- Z


	#Get starting values for the parameters in the model
	optval<-optim(par=c(1, 0.1), fn=loglik.SEBLUP,  gr=gr.SEBLUP,
		direst=direst, X=X, desvar=desvar, 
		Z=Z, W=W, 
		#method="Nelder-Mead"
		method="L-BFGS-B", 
		lower= c(0,-.999),
		upper=c(Inf, .999)#,
		#control=list(factr=.001)
	)


	#Variance of the residuals as initial point     
	sigma2u<-optval$par[1]
	rho<-optval$par[2]

	sigma2urho<-matrix(NA, nrow=2, ncol=1)
	sigma2urho[1,1]<-sigma2u
	sigma2urho[2,1]<-rho

	s<-matrix(NA, nrow=2, ncol=1)
	Isigma2urho<-matrix(NA, nrow=2, ncol=2)


	ZtZ<-Z%*%t(Z)
	WtW<-W%*%t(W)

	niter<-0
	dif<-tol+1

	while((niter<maxiter) & (dif > tol))
	{
		IrhoW<- diag.spam(1, k)-rho*W
		C<-IrhoW%*%t(IrhoW)
		slvC<-solve(C)
		G<-sigma2u*slvC
		GtZ<-G%*%t(Z)
		V<-R + Z%*%GtZ

		Vinv<-as.spam(solve(as.matrix(V))) #FIXME: Use spam only
		VinvX<-Vinv%*%X
		slvXVinvX<-as.spam(solve(as.matrix(t(X)%*%VinvX))) #FIXME

		coeff<-slvXVinvX%*%t(VinvX)%*%direst
		ftd<- X%*%coeff
		rsd<-direst-ftd

		#VinvZtZ<-Vinv%*%ZtZ

		VZCZ<-Vinv%*%Z%*%slvC%*%t(Z)

		s[1,1]<-(-.5)*sum(diag(VZCZ)) +.5*t(rsd)%*%(VZCZ%*%Vinv)%*%rsd

		A<- (-sigma2u)*slvC%*%(2*rho*WtW-2*W)%*%(slvC)
		auxs2<-Vinv%*%Z%*%A%*%t(Z)
		s[2,1]<- (-.5)*sum(diag(auxs2)) +.5*t(rsd)%*%(auxs2%*%Vinv)%*%rsd

		#print(s)

		VZAZ<-Vinv%*%Z%*%A%*%t(Z)

		Isigma2urho[1,1]<- .5*sum(diag(VZCZ%*%VZCZ))
		Isigma2urho[1,2]<- .5*sum(diag(VZCZ%*%VZAZ))
		Isigma2urho[2,1]<- .5*sum(diag(VZAZ%*%VZCZ))
		Isigma2urho[2,2]<- .5*sum(diag(VZAZ%*%VZAZ))

		#print(Isigma2urho)

		sigma2urho<-sigma2urho+solve(Isigma2urho)%*%s

		niter<-niter+1

		dif<-det(Isigma2urho)  # (abs(sigma2u-sigma2urho[1,1]) +abs(rho-sigma2urho[2,1]))/2

		sigma2u<-sigma2urho[1,1]
		rho<-sigma2urho[2,1]

	#print(c(sigma2u, rho))
	}

	#Estimates of the random effects
	randeff<-G%*%t(Z)%*%Vinv%*%(direst-ftd)


	m<-diag.spam(1,k)# m is called b_i^t in Pratesi & Salvati (2008)

	g1<-diag( G-GtZ%*%Vinv%*%t(GtZ) )


	CZVZ<-slvC%*%t(Z)%*%Vinv%*%Z
	AZVZ<-A%*%t(Z)%*%Vinv%*%Z

	g1grad<-matrix(c(
		diag(slvC-(sigma2u*CZVZ%*%slvC-sigma2u*sigma2u*CZVZ%*%CZVZ%*%slvC+sigma2u*CZVZ%*%slvC )) ,
		diag( A-(sigma2u*AZVZ%*%slvC-sigma2u*sigma2u*CZVZ%*%AZVZ%*%slvC+sigma2u*CZVZ%*%A ) ) ),
		nrow=2, byrow=TRUE
		)

	bsigma2urho<- (.5/k)*solve(Isigma2urho)%*%matrix(
			c(sum(diag(slvXVinvX%*%t(X)%*%(-VZCZ%*%VinvX))),
			sum(diag(slvXVinvX%*%t(X)%*%(-VZAZ%*%VinvX)))),
			, ncol=1)

	g1biascorr<-sapply(1:k, function(i){
		g1g<-matrix(c(g1grad[1,i], 
			g1grad[2,i]), ncol=1)

			t(bsigma2urho)%*%g1g

			})

	g2<- sapply(1:k, function(i)
                {
                        tb<-t(m[,i])%*%GtZ%*%VinvX
                        td<-X[i,]-tb

                        (td%*%slvXVinvX%*%t(td))[1,1]
                }) 


	VinvZ<-Vinv%*%Z
	
	g3<-sapply(1:k, function(i)
	{
		delta<-rbind( 
	t(m[,i])%*%( slvC%*%t(VinvZ)+GtZ%*%(-VinvZ%*%slvC%*%t(VinvZ))),
	t(m[,i])%*%( A%*%t(VinvZ) + GtZ%*%(-VinvZ%*%A%*%t(VinvZ)) ) 
		)

		sum(diag(  delta%*%V%*%t(delta)%*%solve(Isigma2urho) ))
	})


	g4<-NULL #FIXME!!! Not computed at the moment


	res<-list(
                coefficients=coeff,
                residuals=as.vector((direst-ftd)[,1]),
                #effects =NULL,
                #rank=NULL,
                fitted.values=as.vector(ftd[,1]),
                randeff=as.vector(randeff[,1]),
                varcoeff=slvXVinvX,
                varsigma2urho=solve(Isigma2urho),
                sigma2u=sigma2u,
		rho=rho,
                g1=g1,
                g2=g2,
                g3=g3,
                g3=g4,
                g1biascorrd=g1biascorr,
                mse= g1-g1biascorr+g2+2*g3,
                eblup=as.vector(ftd[,1]+randeff[,1]),
	loglik=-loglik.SEBLUP(par=c(sigma2u, rho), direst, X, desvar, Z, W),
		method="ML"
        )
        class(res)<-"SEBLUP"

	return(res)

}


SEBLUP.REML<-function(direst, X, desvar, Z=NULL, W, tol=10e-5, maxiter=50)
{
	k<-nrow(direst)

	R<-diag.spam(desvar[,1], k)

	if(is.null(Z))
		Z=diag.spam(1, k)

	l<- t(X)
	m<- Z

	#Get starting values for the parameters in the model
	optval<-optim(par=c(1, 0.1), fn=loglik.SEBLUP,  gr=gr.SEBLUP,
		direst=direst, X=X, desvar=desvar, 
		Z=Z, W=W, 
		#method="Nelder-Mead"
		method="L-BFGS-B", 
		lower= c(0,-.999),
		upper=c(Inf, .999)#,
		#control=list(factr=.001)
	)


	#Variance of the residuals as initial point     
	sigma2u<-optval$par[1]
	rho<-optval$par[2]

	sigma2urho<-matrix(NA, nrow=2, ncol=1)
	sigma2urho[1,1]<-sigma2u
	sigma2urho[2,1]<-rho

	s<-matrix(NA, nrow=2, ncol=1)
	Isigma2urho<-matrix(NA, nrow=2, ncol=2)


	ZtZ<-Z%*%t(Z)
	WtW<-W%*%t(W)

	niter<-0
	dif<-tol+1

	while((niter<maxiter) & (dif > tol))
	{
		IrhoW<- diag.spam(1, k)-rho*W
		C<-IrhoW%*%t(IrhoW)
		slvC<-solve(C)
		G<-sigma2u*slvC
		GtZ<-G%*%t(Z)
		V<-R + Z%*%GtZ

		Vinv<-as.spam(solve(as.matrix(V))) #FIXME: Use spam only
		VinvX<-Vinv%*%X
		slvXVinvX<-as.spam(solve(as.matrix(t(X)%*%VinvX))) #FIXME

		P<-Vinv-VinvX%*%slvXVinvX%*%t(VinvX)

		coeff<-slvXVinvX%*%t(VinvX)%*%direst
		ftd<- X%*%coeff
		rsd<-direst-ftd

		#VinvZtZ<-Vinv%*%ZtZ

		PZCZ<-P%*%Z%*%slvC%*%t(Z)

		s[1,1]<-(-.5)*sum(diag(PZCZ)) +.5*t(direst)%*%(PZCZ%*%P)%*%direst

		A<- (-sigma2u)*slvC%*%(2*rho*WtW-2*W)%*%(slvC)
		auxs2<-P%*%Z%*%A%*%t(Z)
		s[2,1]<- (-.5)*sum(diag(auxs2)) +.5*t(direst)%*%(auxs2%*%P)%*%direst

		#print(s)

		PZAZ<-P%*%Z%*%A%*%t(Z)

		Isigma2urho[1,1]<- .5*sum(diag(PZCZ%*%PZCZ))
		Isigma2urho[1,2]<- .5*sum(diag(PZCZ%*%PZAZ))
		Isigma2urho[2,1]<- .5*sum(diag(PZAZ%*%PZCZ))
		Isigma2urho[2,2]<- .5*sum(diag(PZAZ%*%PZAZ))

		#print(Isigma2urho)

		sigma2urho<-sigma2urho+solve(Isigma2urho)%*%s

		niter<-niter+1

		dif<-det(Isigma2urho)  # (abs(sigma2u-sigma2urho[1,1]) +abs(rho-sigma2urho[2,1]))/2

		sigma2u<-sigma2urho[1,1]
		rho<-sigma2urho[2,1]

	#print(c(sigma2u, rho))
	}

	#Estimates of the random effects
	randeff<-G%*%t(Z)%*%Vinv%*%(direst-ftd)


	m<-diag.spam(1,k)# m is called b_i^t in Pratesi & Salvati (2008)

	g1<-diag( G-GtZ%*%Vinv%*%t(GtZ) )


	CZVZ<-slvC%*%t(Z)%*%Vinv%*%Z
	AZVZ<-A%*%t(Z)%*%Vinv%*%Z

	g1biascorr<-0

	g2<- sapply(1:k, function(i)
                {
                        tb<-t(m[,i])%*%GtZ%*%VinvX
                        td<-X[i,]-tb

                        (td%*%slvXVinvX%*%t(td))[1,1]
                }) 


	VinvZ<-Vinv%*%Z
	
	g3<-sapply(1:k, function(i)
	{
		delta<-rbind( 
	t(m[,i])%*%( slvC%*%t(VinvZ)+GtZ%*%(-VinvZ%*%slvC%*%t(VinvZ))),
	t(m[,i])%*%( A%*%t(VinvZ) + GtZ%*%(-VinvZ%*%A%*%t(VinvZ)) ) 
		)

		sum(diag(  delta%*%V%*%t(delta)%*%solve(Isigma2urho) ))
	})


	g4<-NULL #FIXME!!! Not computed at the moment


	res<-list(
                coefficients=coeff,
                residuals=as.vector((direst-ftd)[,1]),
                #effects =NULL,
                #rank=NULL,
                fitted.values=as.vector(ftd[,1]),
                randeff=as.vector(randeff[,1]),
                varcoeff=slvXVinvX,
                varsigma2urho=solve(Isigma2urho),
                sigma2u=sigma2u,
		rho=rho,
                g1=g1,
                g2=g2,
                g3=g3,
                g3=g4,
                g1biascorrd=g1biascorr,
                mse= g1-g1biascorr+g2+2*g3,
                eblup=as.vector(ftd[,1]+randeff[,1]),
	loglik=-loglik.SEBLUP(par=c(sigma2u, rho), direst, X, desvar, Z, W),
		method="REML"
        )
        class(res)<-"SEBLUP"

	return(res)



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

#FIXME: What likelihood shall we return? Profile, conditional, etc.?
logLik.SEBLUP <- function(object, ..., conditional=TRUE) {

        if(conditional)
{
           mres<-residuals(object) - object$Z%*%ranef(object)
           V<-diag.spam(object$desvar[,1])

        n<-nrow(mres)
        l<- (-n/2)*determinant(V)$modulus -.5*t(mres)%*%solve(V)%*%mres
        }
        else
        {
		l<-object$loglik
        }


        return(l[1,1])
}

