dat <- read.table('Dugong.txt', header = T)
Likelihood <- function(vector){
  Lik = 0
  alpha=vector[1]
  beta = vector[2]
  gamma = vector[3]
  tau.square = vector[4]
  if(((alpha>1)*(beta>1)*(gamma>0 && gamma<1)
       *(tau.square>0))==0){
         return(1e+3)
       }
  for( i in 1:nrow(dat)){
    mu = (alpha-beta*gamma^(dat$Length[i]))
    cat('mu : ', mu , 'prob : ',(1/(sqrt(2*pi*tau.square))*exp(-(dat$Age[i]-mu)^2/(2*tau.square))),
        '-- dnorm : ',dnorm(dat$Age[i],mu,sd=sqrt(tau.square)), '\n')
    Lik = (sum(Lik, log((1/(sqrt(2*pi*tau.square))*exp(-(dat$Age-mu)^2/(2*tau.square))))))
  }
  return(-Lik)
}

initial <- c(1.01, 1.01, 0.0000001, 0.000001)
MLE <- optim(initial,Likelihood, method = 'L-BFGS-B', lower = c(1.01,1.01,0.0000001,0.00001), upper = c(100,100,1,100))

Likelihood(initial)
Likelihood(MLE$par)
  


Likelihood <- function(vector){
  Lik = 0
  alpha=vector[1]
  beta = vector[2]
  gamma = vector[3]
  tau.square = vector[4]
  if(((alpha>1)*(beta>1)*(gamma>0 && gamma<1)
      *(tau.square>0))==0){
    return(1e+9)
  }
  for( i in 1:nrow(dat)){
    mu = (alpha-beta*gamma^(dat$Length[i]))
    cat('mu : ', mu , 'prob : ',(1/(sqrt(2*pi*tau.square))*exp(-(dat$Age[i]-mu)^2/(2*tau.square))),
        '-- dnorm : ',dnorm(dat$Age[i],mu,sd=sqrt(tau.square)), '\n')
    Lik = (prod(Lik, ((1/(sqrt(2*pi*tau.square))*exp(-(dat$Age-mu)^2/(2*tau.square))))))
  }
  return(-Lik)
}
