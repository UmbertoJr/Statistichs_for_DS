samp = rep(NA,1000)
i=0
for(y in Y){
  Z = rnorm(1)*y
  samp[i]= Z
  i = i+1
}

hist(samp, breaks=100, freq = F)
plot(ecdf(samp))