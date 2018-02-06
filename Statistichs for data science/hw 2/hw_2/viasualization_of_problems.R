require(rjags)
require(coda)
require(xtable)

# loading data ------------------------------------------------------------

#obs 1000
load('data_simmulated_for_neal.RData')
load('network_for_K_15_1000.RData')

# muller rios chain
M = 10000

ind <- read.csv(file = "chain winbugs/muller_rois_index.txt", sep = '\t', header = F)
c <- dim(index)[1]
obj <- read.csv(file = "chain winbugs/muller_rois_chain.txt", sep = '\t', header = F)
chain_muller_Rios <- matrix(obj[,2], c(M,c))


# looking to marginal distributions..
l=c(1,16,23,18)
index <- ind[l,1]
s1 <- chain_muller_Rios[,l[1]]
s2 <- chain_muller_Rios[,l[2]]
s3 <- chain_muller_Rios[,l[3]]
s4 <- chain_muller_Rios[,l[4]]

#save(s1,s2,s3,s4,index,ind,l,file = 'presentation/chains_for_muller_plot.RData')

opt1 <- take_the_weight(ind[l[1],1])
opt2 <- take_the_weight(ind[l[2],1])
opt3 <- take_the_weight(ind[l[3],1])
opt4 <- take_the_weight(ind[l[4],1])

par(mfrow=c(2,2), mar=c(2.5,2.5,1.5,2.5))
plot(density(s1), main = ind[l[1],1])
points(opt1,0, col='blue')
plot(density(s2), main = ind[l[2],1])
points(opt2,0, col='blue')
plot(density(s3), main = ind[l[3],1])
points(opt3,0, col='blue')
plot(density(s4), main = ind[l[4],1])
points(opt3,0, col='blue')
opt <- c(opt1,opt2,opt3,opt4)
p_val1 <- sum(s1<opt1)/length(s1)
p_val2 <- sum(s2<opt2)/length(s2)
p_val3 <- sum(s3<opt3)/length(s3)
p_val4 <- sum(s4<opt4)/length(s4)
p_val <- c(min(p_val1,1-p_val1),min(p_val2,1-p_val2),min(p_val3,1-p_val3),min(p_val4,1-p_val4))
table  <- cbind(rbind(HPDinterval(as.mcmc(s1)),HPDinterval(as.mcmc(s2)),
                      HPDinterval(as.mcmc(s3)),HPDinterval(as.mcmc(s4))),opt,p_val)
print(xtable(table))

# Trace plots
par(mfrow=c(2,2), mar=c(2.5,2.5,1.5,2.5))
ts.plot(s1, main = ind[l[1],1])
ts.plot(s2, main = ind[l[2],1])
ts.plot(s3, main = ind[l[3],1])
ts.plot(s4, main = ind[l[4],1])


# Running means
ts.plot(cumsum(s1)/seq_along(s1))
ts.plot(cumsum(sim)/seq_along(s2))
ts.plot(cumsum(sim)/seq_along(sim))
ts.plot(cumsum(sim)/seq_along(sim))

# Auto-correlation functions
acf(sim)

# Cross-Corelation
h <- 5      #1:106
index <- ind[h,1]
index
sim2 <- chain_muller_Rios[,h]
ccf(sim, sim2,50)

# Spectral Analysis with coda Geweke Convergence diagnostic
s1 <- as.mcmc(s1)
Z_score1 <- geweke.diag(s1)
p_v1 <- min(pnorm(Z_score1$z),1 - pnorm(Z_score1$z))

s2 <- as.mcmc(s2)
Z_score2 <- geweke.diag(s2)
p_v2 <- min(pnorm(Z_score2$z),1 - pnorm(Z_score2$z))

s3 <- as.mcmc(s3)
Z_score3 <- geweke.diag(s3)
p_v3 <- min(pnorm(Z_score3$z),1 - pnorm(Z_score3$z))

s4 <- as.mcmc(s4)
Z_score4 <- geweke.diag(s4)
p_v4 <- min(pnorm(Z_score4$z),1 - pnorm(Z_score4$z))

Z_score <- c(Z_score1$z,Z_score2$z,Z_score3$z,Z_score4$z)
p_v <- c(p_v1,p_v2,p_v3,p_v4)
s <- c('Beta1','Beta0','gamma 2,1','gamma 1,2')
table <- cbind(s,Z_score,p_v)
colnames(table) <- c('chain','Z_score','p_value')
row.names(table) <- 1:4
print(xtable(table))

# Heidelberger and Welch Convergence diagnostic
xtable(as.matrix(heidel.diag(s1)))

# posterior prediction distribution
#take l from 107 to 1106
for(obs in 1:10){
  l= 106+ obs
  pred_ind <- ind[l,1]
  sim <- chain_muller_Rios[,l]
  plot(density(sim), main = ind[l,1])
  points(Y[obs],0, col = 'blue')
  p_val <- sum(sim <=Y[obs])/length(sim)
  p_val <- min(p_val, 1- p_val)
  legend('topleft',legend = bquote(p-value ~ "=" ~ .(p_val)),lty = 1:2, cex = 0.5)
  Sys.sleep(0.8)  
}

# visualizing the posterior log-likelihood
l= 1107
pred_ind <- ind[1108:2107,1]
p.log <- rowSums(chain_muller_Rios[,pred_ind])
ts.plot(p.log)




# visualizing distribution for neal ------------------------------------------------
load('neal_chain3_for_15nuerons.RData')
load('neal_chain4_for_15nuerons.RData')
load('neal_chain_huge_chain.RData')

K <- 15
sim <- sim_neal_huge
sim1 <- sim_neal3
sim2 <- sim_neal4

for(k in 1:K){
  plot_and_find_the_mode(sim$Gamma.0,k, obj='Gamma.0')
  Sys.sleep(0.5)
  for(p in 1:P){
    plot_and_find_the_mode(sim$Gamma,k,p, obj='Gamma')
    Sys.sleep(0.5)
  }
}

plot_and_find_the_mode(sim$Beta.0, obj = 'Beta.0')
for(k in 1:K){
  plot_and_find_the_mode(sim$Beta,k, obj = 'Beta')
  Sys.sleep(0.5)
}

for(i in 1:10){
  plot(density(sim$Y.rep[i,,1]))
  points(Y[i],0,col='blue')
  p <- mean(sim$Y.rep[i,,1]< Y[i])
  p_val <- min(p, 1-p)
  print(p_val)
  Sys.sleep(0.7)
}

# Traceplots

ts.plot(sim$Beta[1,,])
ts.plot(sim$Beta.0[1,,])
ts.plot(sim$Gamma.0[2,,])
ts.plot(sim$Gamma[1,1,,])


ts.plot(cbind(sim1$Beta[1,,],sim2$Beta[1,,]))
ts.plot(cbind(sim1$Beta.0[1,,],sim2$Beta.0[1,,]))
ts.plot(cbind(sim1$Gamma.0[2,,],sim2$Gamma.0[2,,]))
ts.plot(cbind(sim1$Gamma[1,1,,],sim2$Gamma[1,1,,]))



# running mean

s1 <- sim$Beta[1,,]
s2 <- sim$Beta.0[1,,]
s3 <- sim$Gamma.0[2,,]
s4 <- sim$Gamma[1,1,,]

ts.plot(cumsum(s1)/seq_along(s1))
ts.plot(cumsum(s2)/seq_along(s2))
ts.plot(cumsum(s3)/seq_along(s3))
ts.plot(cumsum(s4)/seq_along(s4))

# running means two chains

s1_1 <- sim1$Beta[1,,]
s1_2 <- sim1$Beta.0[1,,]
s1_3 <- sim1$Beta[2,,]
s1_4 <- sim1$Gamma[1,1,,]

s2_1 <- sim2$Beta[1,,]
s2_2 <- sim2$Beta.0[1,,]
s2_3 <- sim2$Beta[2,,]
s2_4 <- sim2$Gamma[1,1,,]

ts.plot(cbind(cumsum(s1_1)/seq_along(s1_1),cumsum(s2_1)/seq_along(s2_1)))
ts.plot(cbind(cumsum(s1_2)/seq_along(s1_2),cumsum(s2_2)/seq_along(s2_2)))
ts.plot(cbind(cumsum(s1_3)/seq_along(s1_3),cumsum(s2_3)/seq_along(s2_3)))
ts.plot(cbind(cumsum(s1_4)/seq_along(s1_4),cumsum(s2_4)/seq_along(s2_4)))

# Auto-correlation functions
acf(s1,lag.max = 10000)
acf(s2,10000)
acf(s3,10000)
acf(s4,10000)


# Cross-Corelation
ccf(s1, s2,10000)
ccf(s1, s3,10000)
ccf(s1, s4,10000)

# Spectral Analysis with coda Geweke Convergence diagnostic
s1 <- as.mcmc(s1)
Z_score1 <- geweke.diag(s1)
p_v1 <- min(pnorm(Z_score1$z),1 - pnorm(Z_score1$z))

s2 <- as.mcmc(s2)
Z_score2 <- geweke.diag(s2)
p_v2 <- min(pnorm(Z_score2$z),1 - pnorm(Z_score2$z))

s3 <- as.mcmc(s3)
Z_score3 <- geweke.diag(s3)
p_v3 <- min(pnorm(Z_score3$z),1 - pnorm(Z_score3$z))

s4 <- as.mcmc(s4)
Z_score4 <- geweke.diag(s4)
p_v4 <- min(pnorm(Z_score4$z),1 - pnorm(Z_score4$z))

Z_score <- c(Z_score1$z,Z_score2$z,Z_score3$z,Z_score4$z)
p_v <- c(p_v1,p_v2,p_v3,p_v4)
s <- c('Beta1','Beta0','gamma 2,1','gamma 1,2')
table <- cbind(s,Z_score,p_v)
colnames(table) <- c('chain','Z_score','p_value')
row.names(table) <- 1:4
print(xtable(table,digit = 1))

# Heidelberger and Welch Convergence diagnostic
heidel.diag(s)

# posterior prediction distribution
#take l from 107 to 1106
for(obs in 1:10){
  l= 106+ obs
  pred_ind <- ind[l,1]
  s <- sim$Y.rep[obs,,1]
  plot(density(s), main = ind[l,1])
  points(Y[obs],0, col = 'blue')
  p_val <- sum(s <=Y[obs])/length(s)
  p_val <- min(p_val, 1- p_val)
  legend('topleft',legend = bquote(p-value ~ "=" ~ .(p_val)),lty = 1:2, cex = 0.5)
  Sys.sleep(0.8)  
}

# visualizing the posterior log-likelihood
pred_ind <- ind[1108:2107,1]
p.log <- colSums(sim$p.log[,,1])
ts.plot(p.log)



# third case conjugate pror mixed with flat prior -------------------------
load('chain_conjugate_BNN.RData')
load('chain_conjugate_BNN_2.RData')

sim <- c3$chain
sim1 <- c2
index <- colnames(sim)
# choose l from 1 to 127
l= 12; index[l]
opt <- take_the_weight(index[l]);opt
plot(density(sim[,l]), main = index[l])
points(opt,0, col='blue')
HPDinterval(as.mcmc(sim[,l]))
p_v <- sum(sim[,l]<opt)/length(sim[,l])
print(list('p-value'= min(p_v, 1-p_v)))

#marginal sigma
l= 128; index[l]
plot(density(sim[,l]), main = index[l])
HPDinterval(as.mcmc(sim[,l]))

# Traceplots
s1 <- sim$`Beta 1`
s2 <- sim$`Beta 0`
s3 <- sim$`Gamma 2.1`
s4 <- sim$`Gamma 1.2`


par(mfrow=c(2,2), mar=c(2.5,2.5,1.5,2.5))
ts.plot(s1,main= 'Beta 1')
ts.plot(s2,main= 'Beta 0')
ts.plot(s3,main= 'Gamma.0 2')
ts.plot(s4,main= 'Gamma 1.1')

s1_1 <- s1
s1_2 <- s2
s1_3 <- s3
s1_4 <- s4

s2_1 <- sim1$`Beta 1`
s2_2 <- sim1$`Beta 0`
s2_3 <- sim1$`Gamma 2.1`
s2_4 <- sim1$`Gamma 1.2`

ts.plot(cbind(s1_1,s2_1),main= 'Beta 1')
ts.plot(cbind(s1_2,s2_2),main= 'Beta 0')
ts.plot(cbind(s1_3,s2_3),main= 'Gamma.0 2')
ts.plot(cbind(s1_4,s2_4),main= 'Gamma 1.1')


# running mean

par(mfrow=c(2,2), mar=c(2.5,2.5,1.5,2.5))
ts.plot(cumsum(s1)/seq_along(s1),main= 'Beta 1')
ts.plot(cumsum(s2)/seq_along(s2),main= 'Beta 0')
ts.plot(cumsum(s3)/seq_along(s3),main= 'Gamma.0 2')
ts.plot(cumsum(s4)/seq_along(s4),main= 'Gamma 1.1')

# running means two chains

ts.plot(cbind(cumsum(s1_1)/seq_along(s1_1),cumsum(s2_1)/seq_along(s2_1)))
ts.plot(cbind(cumsum(s1_2)/seq_along(s1_2),cumsum(s2_2)/seq_along(s2_2)))
ts.plot(cbind(cumsum(s1_3)/seq_along(s1_3),cumsum(s2_3)/seq_along(s2_3)))
ts.plot(cbind(cumsum(s1_4)/seq_along(s1_4),cumsum(s2_4)/seq_along(s2_4)))

# Auto-correlation functions
acf(s1,lag.max = 10000)
acf(s2,10000)
acf(s3,10000)
acf(s4,10000)


# Cross-Corelation
ccf(s1, s2,10000)
ccf(s1, s3,10000)
ccf(s1, s4,10000)

# Spectral Analysis with coda Geweke Convergence diagnostic
s <- as.mcmc(s1)
Z_score <- geweke.diag(s);Z_score
min(pnorm(Z_score$z),1 - pnorm(Z_score$z))

# Heidelberger and Welch Convergence diagnostic
heidel.diag(s)

c3$deviance

# functions ---------------------------------------------------------------

take_the_weight <- function(index){
  dim <- length(index)
  for(i in 1:dim){
    i <- as.character(index[i])
    if(length(grep('Beta.0', i))){
      val <- net.1$weights[[1]][[2]][1]
      i <-''
    }
    if(length(grep("Beta", i))){
      v <- as.numeric(gsub(".*([0-9]+).*$", "\\1", index))
      val <- net.1$weights[[1]][[2]][1+v]
      i <-''
    }
    if(length(grep('Gamma.0', i))){
      v <- as.numeric(gsub(".*([0-9]+).*$", "\\1", index))
      val <- net.1$weights[[1]][[1]][1,v]
      i <-''
    }
    if(length(grep("Gamma", i))){
      v <- as.numeric(unlist(strsplit(gsub(".*([0-9]+,[0-9]+).*$", "\\1", index),',')))
      val <- net.1$weights[[1]][[1]][v[2],v[1]]
      i <-''
    }
  }
  return(val)
}

# # Make predictions
# 
# y.hat.mode.2 <- prediction_for_neal_BNN(X,sim = sim)
# y.hat.mean.2 <- prediction_for_neal_BNN_with_mean(X,sim)
# 
# plot(Y,y.hat,xlim = c(-4,5),ylim = c(0,4))
# plot(Y,y.hat.mean.2)
# plot(Y,y.hat.mode.2)
