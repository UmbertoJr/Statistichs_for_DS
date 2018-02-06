
library(truncnorm)

sigma_alpha <- 1000

full.conditional.alpha <- function(start){
  al = (1/(sigma_alpha^2)) + (27/start['tau.square'])
  be = ((start['beta']*(sum(start['gamma']^dat$Age)))   +  
          sum(dat$Length))/ start['tau.square']
  sim_alpha <- 0
  while(sim_alpha<=1){
    sim_alpha <- rtruncnorm(1,mean = be/al , sd = sqrt(1/al), a =1)
  }
  return(sim_alpha)
}


sigma_beta <- 1000

full.conditional.beta <- function(start){
  al = ((sum(start['gamma']^(2*dat$Age))/start['tau.square']) + (1/(sigma_beta^2)))
  be = (start['alpha']*sum(start['gamma']^dat$Age)  -  
          sum(dat$Length * (start['gamma']^dat$Age))) / start['tau.square']
  #cat('mean and variance of beta : ',b/a,1/a,'\n')
  sim_beta <- 0
  while(sim_beta<=1){
    sim_beta <- rtruncnorm(1,mean = be/al , sd = sqrt(1/al), a=1)
  }
  return(sim_beta)
}


gamma.function.ugly.distribution <- function(gam, start){
  alp.n = start['alpha']
  bet.n = start['beta']
  tau.2.n = start['tau.square']
  #cat(alp.n, bet.n, tau.2.n)
  #ret <- (-sum(((bet.n*(gam^dat$Age))*((bet.n*(gam^dat$Age)) - 2*alp.n + 2*dat$Length)))/(2*tau.2.n))
  mu.n = alp.n - bet.n*(gam^dat$Age)
  ret <- exp(- (1/tau.2.n)*0.5*sum((dat$Length - mu.n)^2))
  #cat(ret, '\n')
  return(ret)
}

gamma.function.ugly.distribution <- Vectorize(gamma.function.ugly.distribution)



full.conditional.gamma.1 <- function(start){
  gamma.new <- runif(1)
  check_MH <- runif(1)
  gamm <- start['gamma']
  prob <- min(((gamma.function.ugly.distribution(gamma.new, start))/(gamma.function.ugly.distribution(gamm,start))),1)
  check_MH <- rbinom(1,1,prob = prob)
  if(check_MH){
    gamm <- gamma.new
  }
  return(gamm)
}

full.conditional.gamma.2 <- function(start){
  aa = 1
  bb = 3
  gamma.new <- rbeta(1,shape1 = aa, shape2 = bb)
  gamm <- start['gamma']
  prob <- min(((gamma.function.ugly.distribution(gamma.new, start))/(gamma.function.ugly.distribution(gamm,start)))*dbeta(gamma.new,shape1 = aa,shape2 = bb)/dbeta(gamm,shape1 = aa,shape2 = bb),1)
  check_MH <- rbinom(1,1,prob = prob)
  if(check_MH){
    gamm <- gamma.new
  }
  return(gamm)
}



library(MCMCpack)

al <- 3
be <- 0.5

full.conditional.tau.square <- function(start){
  alpha.n <- start['alpha']
  beta.n <- start['beta']
  gamma.n <- start['gamma']
  a = al + 27/2
  b = ((1/2)*sum((dat$Length - alpha.n + beta.n*(gamma.n^dat$Age))^2)) + be 
  tau.square.sim <- rinvgamma(1,a,b)
  return(tau.square.sim)
}



start<- c('alpha'=1.1, 'beta'=1.01, 'gamma'=0.87,'tau.square'=0.007)

n=10000
simulation1 <- matrix(nrow = n, ncol= 4)
for(i in 1:n){
  start['alpha']= full.conditional.alpha(start)
  start['beta']= full.conditional.beta(start)
  start['gamma'] = full.conditional.gamma.1(start)
  start['tau.square'] = full.conditional.tau.square(start)
  simulation1[i,]=start
}

m <- mcmc(simulation1)
AcceptanceRate(m)

