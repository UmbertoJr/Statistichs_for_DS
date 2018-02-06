library(truncnorm)

# Ages
x = c( 1.0, 1.5, 1.5, 1.5, 2.5, 4.0, 5.0, 5.0, 7.0,
       8.0, 8.5, 9.0, 9.5, 9.5, 10.0, 12.0, 12.0, 13.0,
       13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5)

# Lengths
y = c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
      2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43,
      2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57)


# Hyperparameters

sigma.alpha = 100
sigma.beta = 100
a = 0.001
b = 0.001


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

mh1.rw.unif = function(initialstate, a, loggreekpi, params){
  thetaprop = initialstate+runif(1,min=-a,max=a)
  # acceptance/rejection
  omega = runif(1,min=0,max=1)
  ACCEPT=(omega<min(c(loggreekpi(thetaprop,params)-loggreekpi(initialstate,params),1)))
  theta = initialstate
  if(ACCEPT){
    theta=thetaprop
  }
  return(theta)
}

mh2.beta = function(initialstate, a, loggreekpi, params){
  thetaprop = rbeta(1,shape1 = A,shape2 = B)
  # acceptance/rejection
  omega = runif(1,min=0,max=1)
  ACCEPT=(omega<min(c((loggreekpi(thetaprop,params)-loggreekpi(initialstate,params))*
                      (dbeta(thetaprop,shape1 = A,shape2 = B)/dbeta(initial_state,shape1 = A,shape2 = B)),1)))
  theta = initialstate
  if(ACCEPT){
    theta=thetaprop
  }
  return(theta)
}



n = 10000
alpha.sim = rep(NA,n+1)
beta.sim = rep(NA,n+1)

gamma.sim = rep(NA,n+1)
tau.square.sim = rep(NA,n+1)

alpha.sim[1] = 1
beta.sim[1] = 5
gamma.sim[1] = .8
tau.square.sim[1] = 1




# prova con beta ----------------------------------------------------------

for(g in 1:n) {
  alpha.sim[g+1] = pi.alpha.cond(beta.sim[g],gamma.sim[g],
                                 tau.square.sim[g],x,y,
                                 sigma.alpha)
  beta.sim[g+1] = pi.beta.cond(alpha.sim[g+1],gamma.sim[g],
                               tau.square.sim[g],x,y,
                               sigma.beta)
  params = list(alpha.sim[g+1],beta.sim[g+1],tau.square.sim[g],x,y)
  gamma.sim[g+1] = mh2.beta(gamma.sim[g], .005, prop.log.gamma.cond, params)
  tau.square.sim[g+1] = pi.tau.square.cond(alpha.sim[g+1],beta.sim[g+1],gamma.sim[g+1],x,y,a,b)
  
}
















for(g in 1:n) {
  alpha.sim[g+1] = pi.alpha.cond(beta.sim[g],gamma.sim[g],
                                 tau.square.sim[g],x,y,
                                 sigma.alpha)
  beta.sim[g+1] = pi.beta.cond(alpha.sim[g+1],gamma.sim[g],
                               tau.square.sim[g],x,y,
                               sigma.beta)
  params = list(alpha.sim[g+1],beta.sim[g+1],tau.square.sim[g],x,y)
  gamma.sim[g+1] = mh1.rw.unif(gamma.sim[g], .005, prop.log.gamma.cond, params)
  tau.square.sim[g+1] = pi.tau.square.cond(alpha.sim[g+1],beta.sim[g+1],gamma.sim[g+1],x,y,a,b)
  
}

mean(alpha.sim[2000:n+1])
mean(beta.sim[2000:n+1])
mean(gamma.sim[2000:n+1])
mean(tau.square.sim[2000:n+1])


sim <- cbind(alpha.sim,beta.sim,gamma.sim,tau.square.sim)


m <- mcmc(sim)
AcceptanceRate(m)


B=150
A = ((B-2)*(0.88) - 1)/(0.12)

hist(gamma.sim, freq = 0, probability = 1, breaks = 200, xlim = c(0,1))
curve(dbeta(x,shape1 = A, shape2 = B), add = T, col = 'red')
