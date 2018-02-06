##### Let us start with some MCMC algorithm
##### This is our first Gibbs Sampling algorithm 

#### Let us visualize the full conditionals 
#### for simulating X|Y and Y|X
#### in the case where our target ditribution foy
#### (X,Y) is a (correlated) bivariate gaussian with
#### correlation rho=-0.7

rho=-0.7

library(mvtnorm)

# means vector mu  ...

mu=c(0,0)

# ... and variance-covariance matrix Sigma

Sigma=matrix(c(1,rho,rho,1),ncol=2)

xx=seq(-4,4,0.1)
yy=seq(-4,4,0.1)
temp=cbind(rep(xx,length(yy)),rep(yy,each=length(xx)))
zz=matrix(dmvnorm(temp,mean=mu,sigma=Sigma),nrow=length(xx))

contour(x=xx,y=yy,z=zz) 

sample2=rmvnorm(1000,mean=mu,sigma=Sigma)
points(sample2)


library(rgl)
open3d()
persp3d(x=xx,y=yy,z=zz,xlab="X",ylab="Y",zlab="f(x,y)",main=paste("Bivariate Normal - rho =",rho),col="green3")
grid3d(c("x", "y+", "z"), n =20)



# but also 
persp(xx,yy,zz,zlim=c(0,0.35))
# with titles and a better viewpoint 
persp(xx,yy,zz,zlim=c(0,0.35),theta=30,phi=35,col = "lightblue",xlim=c(-4,4),ylim=c(-4,4))

# we overlap a RED line/band to highlight the initial full-conditionals
# provided that the first initial point at t=0 is x=0.5 y=-0.8

# the full-conditional for x|y 
xx=seq(-4,4,0.1)
yy=c(-0.8,-0.75) 
temp=cbind(rep(xx,2),rep(yy,each=length(xx)))
zz=matrix(dmvnorm(temp,mean=mu,sigma=Sigma),nrow=length(xx))
par(new=T)
persp(xx,yy,zz,zlim=c(0,0.35),theta=30,phi=35,col = "red",xlim=c(-4,4),ylim=c(-4,4))

# the full-conditional for y|x
xx=c(0.4,0.45)
yy=seq(-4,4,0.1)
temp=cbind(rep(xx,each=length(yy)),rep(yy,2))
zz=matrix(dmvnorm(temp,mean=mu,sigma=Sigma),nrow=length(xx),byrow=T)
par(new=T)
persp(xx,yy,zz,zlim=c(0,0.35),theta=30,phi=35,col = "green",xlim=c(-4,4),ylim=c(-4,4))

title(main="Highlighting full-conditionals in the bivariate normal case")


##########################################
#   Gibbs sampling - Bivariate Normal
##########################################

# 

nsim=1000
rho=-0.7

# initialization
theta1.0=-30
theta2.0=30

present.theta1=theta1.0
present.theta2=theta2.0

theta1sim=c(present.theta1)
theta2sim=c(present.theta2)


for(t in 1:nsim){
  
  # generic iteration (updating) at time t 
  
  updated.theta1=rnorm(1,mean=rho*present.theta2,sd=sqrt(1-rho^2))
  updated.theta2=rnorm(1,mean=rho*updated.theta1,sd=sqrt(1-rho^2))
  
  present.theta1=updated.theta1
  present.theta2=updated.theta2
  
  theta1sim=c(theta1sim,present.theta1)
  theta2sim=c(theta2sim,present.theta2)
  
}

thetasim=cbind(theta1sim,theta2sim)

plot(theta1sim,theta2sim);abline(h=0);abline(v=0)

# plot(theta1sim,theta2sim,xlim=c(-5,5),ylim=c(-5,5));abline(h=0);abline(v=0)

# first 10 
# plot(theta1sim[1:10],theta2sim[1:10]);abline(h=0);abline(v=0)
# text(theta1sim[1:10]+1,theta2sim[1:10]+1,label=1:10)

cor(theta1sim,theta2sim)

### now we can write our Gibbs sampling as a function!!

gibbs.biv.norm<-function(iter,theta1,theta2,rho){
  gibbssample<-c()
  for(g in 1:iter){
    theta1<-rnorm(1,mean=rho*theta2,sd=sqrt(1-rho^2))
    theta2<-rnorm(1,mean=rho*theta1,sd=sqrt(1-rho^2))
    gibbssample<-c(gibbssample,theta1,theta2)
  }
  gibbssample<-matrix(gibbssample,nrow=2)
  return(gibbssample)
}



#### NOTE
# to simulate a generic 
# multivariate normal starting from i.i.d multivariate normal
# with mean vector mu ...

mu=c(1,2)

# ...  and variance covariance matrix

Sigma=matrix(c(1,-0.8,-0.8,1),ncol=2)

# using singular value decomposition 
dec=svd(Sigma)

dec
dec$u%*%diag(dec$d)%*%dec$v

# thanks to the matrix A which acts as the "square root of Sigma" or "Sigma^(0.5)"
A=dec$u%*%diag(dec$d)^(0.5)
# so that A%*%t(A) returns Sigma
x=rnorm(2)

x=matrix(rnorm(2000),nrow=2)

y=(A%*%x+mu)

