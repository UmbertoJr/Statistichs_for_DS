
library(truncnorm)


dat$Length <- y
dat$Age <- x


sigma_alpha <- 1000

full.conditional.alpha <- function(start){
  al = (start['tau.square']/(sigma_alpha^2) + 27) /(start['tau.square'])
  be = (sum(dat$Length + start['beta']*start['gamma']^dat$Age))/start['tau.square']
  sim_alpha <- 0
  while(sim_alpha<=1){
    sim_alpha <- rtruncnorm(1,mean = be/al , sd = sqrt(1/al),  a = 1)
  }
  return(sim_alpha)
  }


sigma_beta <- 1000

full.conditional.beta <- function(start){
  al = ((sum(start['gamma']^(2*dat$Age))/start['tau.square']) + (1/(sigma_beta^2)))
  be = (sum((start['alpha']- dat$Length)*(start['gamma']^(dat$Age))))/start['tau.square']
  #cat('mean and variance of beta : ',b/a,1/a,'\n')
  sim_beta <- 0
  while(sim_beta<=1){
    sim_beta <- rtruncnorm(1,mean = be/al , sd = sqrt(1/al), a=1)
  }
  return(sim_beta)
}


gamma.function.ugly.distribution <- function(gam){
  alp.n = start['alpha']
  bet.n = start['beta']
  tau.2.n = start['tau.square']
  #cat(alp.n, bet.n, tau.2.n)
  #ret <- (-sum(((bet.n*(gam^dat$Age))*((bet.n*(gam^dat$Age)) - 2*alp.n + 2*dat$Length)))/(2*tau.2.n))
  mu.n = alp.n - bet.n*gam^dat$Age
  ret <- (- (1/tau.2.n)*0.5*sum((dat$Length - mu.n)^2))
  #cat(ret, '\n')
  return(ret)
}




full.conditional.gamma.1 <- function(start){
  gamma.new <- rbeta(1,shape1 = 0.2, shape2 = 20)
  check_MH <- runif(1)
  gamm <- start['gamma']
  if(check_MH <- exp((gamma.function.ugly.distribution(gamma.new)) 
                     - (gamma.function.ugly.distribution(start['gamma'])))){
    gamm <- gamma.new
  }
  return(gamm)
}

full.conditional.gamma.2 <- function(start){
  aa = 0.5
  bb = 330
  gamma.new <- rbeta(1,shape1 = aa, shape2 = bb)
  check_MH <- runif(1)
  gamm <- start['gamma']
  if(check_MH <- ((gamma.function.ugly.distribution(gamma.new)) 
                     - (gamma.function.ugly.distribution(start['gamma'])))
     *dbeta(gamma.new,shape1 = aa,shape2 = bb)/dbeta(start['gamma'],shape1 = aa,shape2 = bb)){
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
  tau.square.sim <- rgamma(1,a,b)
  return(1/tau.square.sim)
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

old.par <- par(mfrow=c(2,2), mar=c(1,1,1,1))

par(old.par)


plot(simulation[,1], type = 'l')
hist(simulation[,1], breaks = 150)

plot(simulation[,2], type = 'l')
hist(simulation[,2], breaks = 150)

library(mcmc)

plot(simulation[,3], type = 'l')
hist(simulation[,3], breaks = 150, probability = 1)
curve(dbeta(x,shape1 =0.5 ,shape2 = 330), add = T, col= 'red')



plot(simulation[,4], type = 'l')
hist(simulation[,4], breaks = 150)

summary(simulation[,1])

summary(simulation[,2])

summary(simulation[,3])

summary(simulation[,4])

plot(cumsum(simulation[,1])/seq_along(simulation[,1]), type='l')
  
plot(cumsum(simulation[,2])/seq_along(simulation[,2]), type='l')


plot(cumsum(simulation[,3])/seq_along(simulation[,3]), type='l')

plot(cumsum(simulation[,4])/seq_along(simulation[,4]), type='l')



M <- mcmc(data = simulation)

AcceptanceRate(M)
