

library(SAE)
spam.options(eps=.0000000001)

data(seblup)

direst<-ydir1
X<-Xpop
desvar<-vardir

d<-data.frame(direst=direst, covariate=Xpop[,2], desvar=desvar)





eblupml<-EBLUP.ML(direst, X, desvar)


eblupml<-EBLUP(direst~covariate, ~desvar, d)



eblupreml<-EBLUP(direst~covariate, ~desvar, d, method="REML")


summary(eblupreml$randeff - eblupml$randeff)

plot(eblupreml$randeff , eblupml$randeff)
abline(0,1)


plot(eblupml$mse, eblupreml$mse)
abline(0,1)


seblupml<-SEBLUP.ML(direst, X, desvar, NULL, W)


seblupreml<-SEBLUP.REML(direst, X, desvar, NULL, W)



#Test unit level model
library(SAE)
spam.options(eps=.0000000001)

data(seblup)

source("/home/vrubio/Imperial/ESRC-Project/Code/RSAE/SAE2/R/SEBLUP_unit.R")

d<-data.frame(y=S[,2], covariate=XS[,1], municip=XS[,2]) 
dpop<-data.frame(covariate=Xpop[,2], municip=1:42, y=as.numeric(by(d$y, d$municip, mean)) )


#Make a bigger data set
d2<-d[rep(1:nrow(d), times=100),]

spam.options(eps=.000001)
unitseblup<-seblup(y~covariate, ~municip,  data=d, datapop=dpop, W=W)
#unitseblup<-seblup(y~covariate, ~municip,  data=d2, datapop=dpop, W=W, tol=spam.options()$eps )

plot(dpop$y, unitseblup$eblup)
abline(0,1)


actmeans<-as.numeric(by(y[,2], y[,1], mean))

plot(actmeans, unitseblup$eblup, pch=19)
points(actmeans, seblupreml$eblup)
points(actmeans, eblupml$eblup, pch=3)
legend(0, 120, bty="n", c("USEBLUP", "SEBLUP", "EBLUP"), pch=c(19, 1, 3))
abline(0,1)

print(sum( (actmeans-unitseblup$eblup)^2)/284)
print(sum( (actmeans-seblupreml$eblup)^2)/284)
print(sum( (actmeans-eblupml$eblup)^2)/284)

