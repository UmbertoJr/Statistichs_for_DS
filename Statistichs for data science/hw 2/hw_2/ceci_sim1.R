library(truncnorm)

pi.alpha.cond = function(beta,gamma,tau.square,x,y,sigma.alpha) {
  d = tau.square+(length(y)*sigma.alpha^2)
  mu = (sum(y+beta*gamma^x)*sigma.alpha^2)/d
  sigma = sqrt((tau.square*sigma.alpha^2)/d)
  alpha = rtruncnorm(1,1,Inf,mu,sigma)
  return(alpha)
}

pi.beta.cond = function(alpha,gamma,tau.square,x,y,sigma.beta) {
  d = tau.square+(sum(gamma^(2*x))*sigma.beta^2)
  mu = (sum(alpha*(gamma^x)-(y*(gamma^x)))*sigma.beta^2)/d
  sigma = sqrt((tau.square*sigma.beta^2)/d)
  beta = rtruncnorm(1,1,Inf,mu,sigma)
  return(beta)
}

pi.tau.square.cond = function(alpha,beta,gamma,x,y,a,b) {
  a.new = (length(y)/2)+a
  b.new = (2*b+sum((y-alpha+(beta*gamma^x))^2))/2
  tau.square = rinvgamma(1,a.new,b.new)
  return(tau.square)
}

prop.log.gamma.cond = function(gamma,param) {
  alpha = param[[1]]
  beta = param[[2]]
  tau.square = param[[3]]
  x = param[[4]]
  y = param[[5]]
  p = (-1/(2*tau.square))*sum((beta^2*gamma^(2*x))+(2*y*beta*gamma^x)-(2*alpha*beta*gamma^x))
  if (is.nan(p)) p = -Inf
  return(p)
}


mh = function(initialstate, a, loggreekpi, params){
  thetaprop = initialstate+runif(1,min=-a,max=a)
  # acceptance/rejection
  omega = runif(1,min=0,max=1)
  ACCEPT=(omega<min(c(exp(loggreekpi(thetaprop,params)-loggreekpi(initialstate,params)),1)))
  theta = initialstate
  if(ACCEPT){
    theta=thetaprop
  }
  return(theta)
}



sigma.alpha = 100
sigma.beta = 100
a = 0.001
b = 0.001

n = 10000


alpha.sim = rep(NA,n+1)
beta.sim = rep(NA,n+1)
gamma.sim = rep(NA,n+1)
tau.square.sim = rep(NA,n+1)

alpha.sim[1] = 1
beta.sim[1] = 5
gamma.sim[1] = .8
tau.square.sim[1] = 1

for(g in 1:n) {
  alpha.sim[g+1] = pi.alpha.cond(beta.sim[g],gamma.sim[g],
                                 tau.square.sim[g],x,y,
                                 sigma.alpha)
  beta.sim[g+1] = pi.beta.cond(alpha.sim[g+1],gamma.sim[g],
                               tau.square.sim[g],x,y,
                               sigma.beta)
  params = list(alpha.sim[g+1],beta.sim[g+1],tau.square.sim[g],x,y)
  gamma.sim[g+1] = mh(gamma.sim[g], .005, prop.log.gamma.cond, params)
  tau.square.sim[g+1] = pi.tau.square.cond(alpha.sim[g+1],beta.sim[g+1],gamma.sim[g+1],x,y,a,b)
  
}


sim <- cbind(alpha.sim,beta.sim,gamma.sim,tau.square.sim)


m <- mcmc(sim)
AcceptanceRate(m)




mean(alpha.sim[2000:n+1])
mean(beta.sim[2000:n+1])
mean(gamma.sim[2000:n+1])
mean(tau.square.sim[2000:n+1])