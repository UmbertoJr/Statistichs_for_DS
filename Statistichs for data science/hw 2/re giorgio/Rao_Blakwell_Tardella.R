## R code for the example of Rao-Blackwell-ization of 
## a monte Carlo Estimate

g=5
mcsize=10000

T=rt(mcsize,df=g)
I_hat = mean(exp(-T^2))

Y = rgamma(mcsize,rate=g/2,shape=g/2)
cond.expectation=1/sqrt(1+2/Y)
I_RB_hat= mean(cond.expectation)



I_hat_sim=rep(NA,replicate)

replicate=500
I_hat_sim=rep(NA,replicate)
I_RB_hat_sim=rep(NA,replicate)

for(i in 1:replicate){
  
  T=rt(mcsize,df=g)
  I_hat_sim[i]=mean(exp(-T^2))
  
  Y = rgamma(mcsize,rate=g/2,shape=g/2)
  cond.expectation=1/sqrt(1+2/Y)
  I_RB_hat_sim[i]=mean(cond.expectation) 
  
}

plot(density(I_RB_hat_sim))
lines(density(I_hat_sim))


