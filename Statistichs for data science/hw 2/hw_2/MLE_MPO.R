dat <- read.table('Dugong.txt', header = T)



res <- data.frame()
i=1
for(a in 1:3){
  for(b in 1:3){
    for(g in 1:3){
      for(t in 1:3){
        i=i+1
        cat(i,' ',Likelihood(c(1.1+(a-1)*100, 1.1+(b-1)*90, 0.2+(g-1)*0.2, 3+(t-1)*100)),'\n')
        cat(c(1.1+(a-1)*100, 1.1+(b-1)*100, 0.2+(g-1)*0.2, 3+(t-1)*100),'\n')
        MLE <- optim(c(1.1+(a-1)*100, 1.1+(b-1)*100, 0.2+(g-1)*0.2, 3+(t-1)*100),Likelihood)
        res <- rbind(res,c(MLE$par,Likelihood(MLE$par)))
      }
    }
  }
}

res


opt = res[res[,6]==min(res[,6])][2:6]
cat('MLE search alpha =',opt[1],'\n',
    'MLE search beta =',opt[2],'\n',
    'MLE search gamma =',opt[3],'\n',
    'MLE search tau.square =',opt[4],'\n',
    '- log(Likelihood) optim =',Likelihood(opt[1:5]),'\n')




Prior <- function(vector){
  alpha=vector[1]
  beta = vector[2]
  gamma = vector[3]
  tau.square = vector[4]
  alpha.p <- dnorm(alpha,0,1e+10)*(alpha>=1)
  beta.p <- dnorm(beta,0,1e+10)*(beta>=1)
  gamma.p <- dunif(gamma,0,1)
  tau.square <- dinvgamma(tau.square,1,1)*(tau.square>=0)
  join <- log(alpha.p)+log(beta.p)+log(gamma.p)+log(tau.square)
  return(-join)
}

Posterior <- function(vector){Prior(vector) + Likelihood_3(vector)}

#Posterior(c(1.500042819, 1.552780102, 0.003414928, 0.001686139))
#Posterior(c(1.500042821, 1.000000018, 0.003414923, 0.001686130))
#MLE <- optim(c(1.1,1.14,0.26,0.0005),Posterior, control = list(fnscale=-1))
MPO <- optim(c(1.1, 1.1, 0.2, 3),Posterior)

cat('Maximum a Posterior alpha.hat =',MLE$par[1],'\n',
    'Maximum a Posterior beta.hat =',MLE$par[2],'\n',
    'Maximum a Posterior gamma =',MLE$par[3],'\n',
    'Maximum a Posterior tau.square =',MLE$par[4],'\n',
    '-log(Posterior) local optim =',Posterior(MPO$par),'\n',
    '-log(Posterior) initial value =',Posterior(c(1.1, 1.1, 0.2, 3)))

mu.func <- function(x, vector= vector){
  alpha=vector[1]
  beta = vector[2]
  gamma = vector[3]
  tau.square = vector[4]
  mu = (alpha-beta* (gamma^(x)))
  return(mu)
}

plot(dat$Age,dat$Length, col='black',xlab='Age',ylab='Length')
curve(mu.func(x, vector = opt[1:4]), add = TRUE, col='blue')
curve(mu.func(x, vector = opt_mpo[1:4]), col = 'red', add= TRUE)



i=1
MPO <- optim(c(1.1, 1.1, 0.2, 3),Posterior)
res_mpo <- matrix(c(i,MPO$par,Posterior(MPO$par)), ncol = 6)

for(a in 1:3){
  for(b in 1:3){
    for(g in 1:3){
      for(t in 1:3){
        i=i+1
        cat(i,' ',Posterior(c(1.1+(a-1)*1, 1.1+(b-1)*1, 0.2+(g-1)*0.1, 3+(t-1)*1)),'\n')
        cat(c(1.1+(a-1)*1, 1.1+(b-1)*1, 0.2+(g-1)*0.1, 3+(t-1)*1),'\n')
        MPO <- optim(c(1.1+(a-1)*1, 1.1+(b-1)*1, 0.2+(g-1)*0.1, 3+(t-1)*1),Posterior)
        res_mpo <- rbind(res_mpo,c(i, MPO$par, Posterior(MPO$par)))
      }
    }
  }
}

save(res_mpo, file = 'risultati_mpo.RData')
  

opt_mpo = unlist(res_mpo[res_mpo[,6]==min(res_mpo[,6]),][2:6])
cat('MPO search alpha =',opt_mpo[1],'\n',
    'MPO search beta =',opt_mpo[2],'\n',
    'MPO search gamma =',opt_mpo[3],'\n',
    'MLE search tau.square =',opt_mpo[4],'\n',
    '- log(Posterior) optim =',Posterior(opt_mpo[1:4]),'\n')























Likelihood_1 <- function(vector){
  Lik = 0
  alpha=vector[1]
  beta = vector[2]
  gamma = vector[3]
  tau.square = vector[4]
  if(((alpha>1)*(beta>1)*(gamma>0 && gamma<1)
                    *(tau.square>0))==0){
    Lik = 1e+10
    return(Lik)
  }
  for( i in 1:nrow(dat)){
    mu = (alpha-beta*gamma^(dat$Length[i]))
    l = log(tau.square^(-0.5)) - ((dat$Age[i] - mu)^2 / (2*tau.square))
    Lik = (sum(Lik, l))
  }
  return(-Lik)
}

MLE <- optim(c(1.1, 1.1, 0.2, 3),Likelihood_1)


Likelihood_3 <- function(vector){
  Lik = 0
  alpha=vector[1]
  beta = vector[2]
  gamma = vector[3]
  tau.square = vector[4]
  if(((alpha>1)*(beta>1)*(gamma>0 && gamma<1)
      *(tau.square>0))==0){
    Lik = 1e+10
    return(Lik)
  }
  for( i in 1:nrow(dat)){
    mu = (alpha-beta* (gamma^(dat$Age[i])))
    if(mu <0){
      Lik = 1e+10
      return(Lik)
    }
    l = log(tau.square^(-0.5)) - ((dat$Length[i] - mu)^2 / (2*tau.square))
    Lik = (sum(Lik, l))
  }
  return(-Lik)
}

MLE <- optim(c(1.1, 1.1, 0.2, 3),Likelihood_3)


i=1
MLE <- optim(c(1.1, 1.1, 0.2, 3),Likelihood_3)
res <- matrix(c(i,MLE$par,Likelihood_3(MLE$par)), ncol = 6)

for(a in 1:3){
  for(b in 1:3){
    for(g in 1:3){
      for(t in 1:3){
        i=i+1
        cat(i,' ',Likelihood_3(c(1.1+(a-1)*1, 1.1+(b-1)*1, 0.2+(g-1)*0.1, 3+(t-1)*1)),'\n')
        cat(c(1.1+(a-1)*1, 1.1+(b-1)*1, 0.2+(g-1)*0.1, 3+(t-1)*1),'\n')
        MLE <- optim(c(1.1+(a-1)*1, 1.1+(b-1)*1, 0.2+(g-1)*0.1, 3+(t-1)*1),Likelihood_3)
        res <- rbind(res,c(i, MLE$par, Likelihood_3(MLE$par)))
      }
    }
  }
}

res

opt = unlist(res[res[,6]==min(res[,6]),][2:6])
cat('MLE search alpha =',opt[1],'\n',
    'MLE search beta =',opt[2],'\n',
    'MLE search gamma =',opt[3],'\n',
    'MLE search tau.square =',opt[4],'\n',
    '- log(Likelihood) optim =',Likelihood_3(opt[1:4]),'\n')






mu_find <- function(vector){
  alpha=vector[1]
  beta = vector[2]
  gamma = vector[3]
  tau.square = vector[4]
  mu = (alpha-beta* pow(gamma,(dat$Age)))
  return(mu)
}


Likelihood_2(c(1.1,1.1,0.1,3))
Likelihood_1(c(1.1,1.1,0.1,3))

Likelihood_2 <- function(vector){
  Lik = 0
  alpha=vector[1]
  beta = vector[2]
  gamma = vector[3]
  tau.square = vector[4]
  if(((alpha>1)*(beta>1)*(gamma>0 && gamma<1)
      *(tau.square>0))==0){
    Lik = 1e+10
    return(Lik)
  }
  
  mu = (alpha-beta*gamma^(dat$Length))
  l = log(tau.square^(-0.5)) - (sum((dat$Age - mu)^2) / (2*tau.square))
  Lik = (sum(Lik, l))
  return(-Lik)
}

