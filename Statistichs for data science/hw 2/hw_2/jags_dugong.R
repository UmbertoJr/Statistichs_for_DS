### Bayesian linear mixed effects model (random intercepts + slopes model)

setwd("~/Dropbox/2017-Stat4DS2-CS/Lab/LMM")

# Simulated longitudinal data

n=75	# number of observed individuals
M=5		# number of measurements 
N=n*M	# balanced design (M is constant over the individuals)

occasions=c(0,2,4,8,10)
set.seed(849367)
gender=sample(c("F","M"),size=n,replace=TRUE,prob=c(0.6,0.4))

ID=rep(1:n,each=M)		# grouping variable
times.long=rep(occasions,n)  # time covariate
gender.long=rep(gender,each=M)	# gender covariate

X=model.matrix(~1+times.long+gender.long)	# design matrix of fixed effects
Z=X[,1:2]	# random intercept + random slopes (random effect of the time)

# Set the fixed effects parameters...

beta0=40
beta1=12
beta2=25

beta=c(beta0,beta1,beta2)
K=length(beta)

# ...and the variance parameters

sigma.e=10

sigma.b0=8
sigma.b1=6

# Two-step procedure to simulate the reponses:

# 1) Generate the random effects from their marginal distributions

set.seed(849367)
b0=rnorm(n=n,mean=0,sd=sigma.b0)

set.seed(849367)
b1=rnorm(n=n,mean=0,sd=sigma.b1)

b=cbind(b0,b1)[ID,]

# 2) Generate the outcomes from their conditional distributions

set.seed(849367)
y=rnorm(n=N,mean=X%*%beta+rowSums(Z*b),sd=sigma.e)

# Plot of the individual trajectories alltogether

matplot(matrix(y,nrow=M,ncol=n),type="l",ylab="y",xlab="Occasion t",ylim=c(-50,300),
        main="Individual trajectories")

matplot(matrix(y,nrow=M,ncol=n),type="l",ylab="y",xlab="Occasion t",ylim=c(-50,300),
        main="Individual trajectories",col=c("violet","blue")[(gender=="M")+1])

# Plot of the individual trajectories

library("lattice")

xyplot(y~times.long|factor(ID),type="l")

# Let us fit the random intercept + slope model

# MLE analysis using R

library(lme4)

lme.REML=lmer(y~times.long+gender.long+(1|ID)+(0+times.long|ID))
fixef(lme.REML)
ranef(lme.REML)
coef(lme.REML)
print(lme.REML)

par(mfrow=c(1,2))

plot(b0,ranef(lme.REML)$ID[,1],xlab=expression(b[0][i]),ylab=expression(hat(b)[0][i]),
     cex.lab=.8)
abline(a=0,b=1)
plot(b1,ranef(lme.REML)$ID[,2],xlab=expression(b[1][i]),ylab=expression(hat(b)[1][i]),
     cex.lab=.8)
abline(a=0,b=1)

# Bayesian analysis using JAGS

# Model specification

cat("model {
    
    # Priors
    
    for(k in 1:K){
    beta[k]~dnorm(0,0.001)
    }	
    
    tau.e~dgamma(0.001,0.001)
    sigma.e <- 1/sqrt(tau.e)
    
    tau.b0~dgamma(0.001,0.001)
    sigma.b0 <- 1/sqrt(tau.b0)
    
    tau.b1~dgamma(0.001,0.001)
    sigma.b1 <-1 /sqrt(tau.b1)
    
    # Statistical (conditional) model
    
    for(i in 1:n){
    b0[i]~dnorm(0,tau.b0)
    b1[i]~dnorm(0,tau.b1)
    for(t in 1:M){
    mu[i,t] <- beta[1] + b0[i] + (beta[2] + b1[i])*times[i,t] + beta[3]*gender[i,t]
    y[i,t]~dnorm(mu[i,t],tau.e)
    }
    }
    
    }",file="lmm-true.model.txt",fill=TRUE)

# Data preparation

data.input=list(n=n,M=M,K=K,
                y=matrix(y,nrow=n,ncol=M,byrow=TRUE),
                times=matrix(times.long,nrow=n,ncol=M,byrow=TRUE),
                gender=matrix(X[,"gender.longM"],nrow=n,ncol=M,byrow=TRUE)) 

# Initialization values

inits=list(
  list(b0=ranef(lme.REML)$ID[,1],b1=ranef(lme.REML)$ID[,2],
       beta=as.numeric(fixef(lme.REML)),
       tau.e=1/sigma(lme.REML)^2,
       tau.b0=1/as.data.frame(VarCorr(lme.REML))[1,"sdcor"]^2,
       tau.b1=1/as.data.frame(VarCorr(lme.REML))[2,"sdcor"]^2)
)

# Parameters to monitor

params=c("beta","sigma.e","sigma.b0","sigma.b1")

# Let us call JAGS

library(R2jags)

true.model.jags=jags(data=data.input,inits=inits,
                     parameters.to.save=params,
                     model.file="lmm-true.model.txt",
                     DIC=TRUE,n.chains=1,n.iter=11000,n.burnin=1000,n.thin=1)
print(true.model.jags,digits=3)
traceplot(true.model.jags)
print(true.model.jags$BUGSoutput$mean,digits=3)

cat("model {
    
    # Priors
    
    for(k in 1:K){
    beta[k]~dnorm(0,0.001)
    }	
    
    tau.e~dgamma(0.001,0.001)
    sigma.e <- 1/sqrt(tau.e)
    
    tau.b0~dgamma(0.001,0.001)
    sigma.b0 <- 1/sqrt(tau.b0)
    
    # Statistical (conditional) model
    
    for(i in 1:n){
    b0[i]~dnorm(0,tau.b0)
    for(t in 1:M){
    mu[i,t] <- beta[1] + b0[i] + beta[2]*times[i,t] + beta[3]*gender[i,t]
    y[i,t]~dnorm(mu[i,t],tau.e)
    }
    }
    
    }",file="lmm-model2.txt",fill=TRUE)

lme.REML2=lmer(y~times.long+gender.long+(1|ID))

# Initialization values

inits2=list(
  list(b0=ranef(lme.REML)$ID[,1],
       beta=as.numeric(fixef(lme.REML)),
       tau.e=1/sigma(lme.REML)^2,
       tau.b0=1/as.data.frame(VarCorr(lme.REML))[1,"sdcor"]^2)
)

# Parameters to monitor

params2=c("beta","sigma.e","sigma.b0")

model2.jags=jags(data=data.input,inits=inits2,
                 parameters.to.save=params2,
                 model.file="lmm-model2.txt",
                 DIC=TRUE,n.chains=1,n.iter=11000,n.burnin=1000,n.thin=1)
print(model2.jags,digits=3)
traceplot(model2.jags)
print(model2.jags$BUGSoutput$mean,digits=3)




