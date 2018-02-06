# last update: 04.20.2017
# (update n.obs after reading real data set!) 

# from Carlin, Gelfand e Smith (1992) 
# Applied Statistics, 41, 389-405
# Hierarchical Bayesian Analysis of Change-Point Problems

# BEST PRACTICE:
# before implementing a GS to perform your Bayesian real data analysis
# with a "new" model
# 
# do first simulate synthetic data from the model (say with n.obs=100) 
# with known TRUE parameter values
# and make inference on the simulated data 

# INGREDIENTS:

# TRUE (fixed, known) parameters

set.seed(1234)
n.obs=100
m.true=30
lambda.true=1.5
phi.true=3.5

# synthetic data

# from a *KNOWN* MODEL 

y.firstperiod=rpois(n=m.true,lambda=lambda.true)
y.secondperiod=rpois(n=n.obs-m.true,lambda=phi.true)

y=c(y.firstperiod,y.secondperiod)

plot(y,xlab="t",type="h")
title=(main="Simulated data")

points(1:m.true,y.firstperiod,pch="x",col="red",font=2)
points((m.true+1):n.obs,y.secondperiod,pch="+",col="blue",font=2)

# NOTICE WE CAN color THIS PLOT because WE HAVE SIMULATED OUR DATA!!

# We would like to approximate (among other things)
# (posterior) expectations of 

# lambda

# phi 

# m



# full conditionals 

# CAREFUL! parameterization of the gamma density in R 

# lambda|y,phi,m ~ Gamma(ip.alpha.rate+m,ip.beta.shape+sum(y))
# phi|y,m,lambda ~ Gamma(ip.a.rate+n-m,ip.b.shape+sum(y))
# m|y,lambda,phi ~ exp(log(lambda)*(newbeta.shape+sum(y[1:m])-1) -lambda*(newalpha.rate+m)
#                      + log(phi)*(ip.b.shape+sum(y[(m+1):n.obs])-1) - phi*(b+n-m)         )  m=1,...,n-1

# fix some (prior) values on the hyperparameters (parameter of the prior distribution)

ip.alpha.rate=0.1
ip.beta.shape=0.2
ip.a.rate=0.1
ip.b.shape=0.7

# what uncertainty on the parameters?


# fix the length of the simulated chains
# and then 
# initialize the corresponding vectors

n.iter=500

lambda=rep(NA,n.iter+1)
phi=rep(NA,n.iter+1)
m=rep(NA,n.iter+1)

m.support=seq(1,n.obs-1)

# initialize the chain at time t=0 [1]

lambda.start=lambda[1]=1.5
phistart=phi[1]=3.5
m.start=m[1]=30

# is it fair? maybe ... well, of course not ... but this is certainly a good starting point!

# otherwise 

lambda.start=lambda[1]=5
phistart=phi[1]=7
m.start=m[1]=60


stat.y.firstperiod=cumsum(y)[-length(y)]
stat.y.secondperiod=sum(y)-stat.y.firstperiod
# for all m in m.support

for(t in 1:n.iter){
  
  # UPDATE POISSON PARAMETER FOR THE FIRST PERIOD 
  
  lambda[t+1]=rgamma(1,shape=ip.beta.shape+sum(y[1:m[t]]),rate=ip.alpha.rate+m[t])
  # lambda[t+1]=rgamma(1,shape=ip.beta.shape+stat.y.firstperiod[m[t]],rate=ip.alpha.rate+m[t])
  
  # UPDATE POISSON PARAMETER FOR THE SECOND PERIOD 
  
  
  phi[t+1]=rgamma(1,shape=ip.b.shape+sum(y[(m[t]+1):n.obs]),rate=ip.a.rate+n.obs-m[t])
  # phi[t+1]=rgamma(1,shape=ip.b.shape+stat.y.secondperiod[m[t]],rate=ip.a.rate+n.obs-m[t])
  
  # UPDATE changepoint 
  # using the full-conditional (discrete) masses
  
  # CAREFUL!! we need to rescale on a proper range in order to avoid overflows
  # related to exponentials
  # (use log-scale!)
  
  logci=log(lambda[t+1])*(stat.y.firstperiod)- lambda[t+1]*(m.support)+ log(phi[t+1])*(stat.y.secondperiod)- phi[t+1]*(n.obs-m.support)
  
  # back to the natural scale
  m.full.conditional.nn=exp(logci-max(logci))
  
  # NOTE that sample(...) does not need probability masses but *positive* masses (up to proportionality constant)
  
  m[t+1]=sample(x=m.support,size=1,prob=m.full.conditional.nn)
  
  
}

jointchain=data.frame(lambda,phi,m)


# FIRST OF ALL: TRACE PLOTS to get a feeling of the influence of the starting point and possible impact of autocorrelation 

plot(lambda,type="l")
acf(lambda)

plot(phi,type="l")
acf(phi)

plot(m,type="l")
acf(m)

# THE (MCMC) SUMMARIES

plot(table(m[]))

# plot(table(m[25000:50001]))


### Now let us explore the result of a simulation 
### from a Markov Chain with stationary distribution 
### corresponding to the (joint) posterior distribution 
### of the parameter of our change-point statistical model of interest

par(mfrow=c(2,2))
hist(m)
plot(lambda,phi)
plot(table(m),type="h")


# NOW LET US PROVIDE SOME INFERENTIAL FINDINGS WITH OUR BAYESIAN ANALYSIS AND COMPARE THEM WITH THE TRUE MODEL PARAMETERS!!


# one way of summarizing the marginal posterior of m|data
prop.table(table(m))

names(which.max(prop.table(table(m))))

# alternatively
mean(m)



# TWO IMPORTANT THINGS 

# Don't forget to (possibly) remove some of the initial T_0 simulations

# Don't forget to assess posterior uncertainty 


# WHAT ABOUT THE OTHER TWO PARAMETERS?


summary(lambda)
mean(lambda)
quantile(lambda,prob=c(0.025,0.975))

quantile(lambda,seq(0.1,0.9,0.1))

# SIMILARLY FOR phi

summary(phi)
mean(phi)
quantile(phi,prob=c(0.025,0.975))

quantile(phi,seq(0.1,0.9,0.1))

# OTHER FE

# axis(1,at=seq(2,100,2),cex=0.5,las=2)
hist(lambda/phi)

mean(lambda/phi)

cor(lambda,phi)



# what do you expect if I change (increase) n.obs ?

# what kind of posterior summaries would you report?






# now *you* can use the real data to carry on your first approximate MCMC Bayesian inference!!!

disasters=read.csv("coal-mining-disasters.csv")
str(disasters)


plot(disasters$year,disasters$n.disasters,main="Number of fatal accidents in UK coal mining sites\n[real data]",xlab="year",ylab="number of accidents",type="h")

# recycle your previous code but .... careful!!

y=disasters$n.disasters
n.obs=length(y)

m.support=seq(1,n.obs-1)
plot(y,type="h")

set.seed(1234)

# re-use the previous code ......