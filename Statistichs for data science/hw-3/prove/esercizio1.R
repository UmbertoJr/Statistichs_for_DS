library(mnormt)

# Simulation size
n= 1e4

# This (M x 7) matrix will store the samples from the joint
sim.samp <- matrix(NA, nrow = M, ncol = 8)
colnames(sim.samp) <- c("sigma1","sigma2","tetha","x1","x2","y1","y2", "DAG")

#Main loop
varcov<- matrix(c(1,0,0,1),2,2)
for (i in 1:n){
  sim.samp[i,1:2]<- rmnorm(1,c(2,7), varcov)
  sim.samp[i,3]<- abs(rnorm(1))
  sim.samp[i,4] <- rnorm(1, mean = sim.samp[i,1], sd = sim.samp[i,3])
  sim.samp[i,5] <- rnorm(1, mean = sim.samp[i,2], sd = sim.samp[i,3])
  sim.samp[i,6] <- rnorm(1, mean = sim.samp[i,1], sd = sim.samp[i,3])
  sim.samp[i,7] <- rnorm(1, mean = sim.samp[i,2], sd = sim.samp[i,3])
  sim.samp[i,8] <- prod(sim.samp[i,1:7])
}


var(sim.samp[,5])
plot(ecdf(sim.samp[,8]),main='ecdf')
abs(-2)
save(sim.samp, file = 'sim_samp.RData')
hist(sim.samp[,8],main = 'DAG distro',breaks = 200)


curve(dnorm(x,mean = mean(sim.samp[,8]) , sd = sd(sim.samp[,8])), add = T, col = "purple", lwd = 4)
      