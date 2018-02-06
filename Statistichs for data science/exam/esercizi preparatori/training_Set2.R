require(boot, quietly = T)

# True interquantila range
iqr.true <- (qt(0.75, df=3) - qt(0.25, df=3))/1.34

# Empirical interquantile range
iqr <- function(x, idx) (quantile(x[idx], 0.75) - quantile(x[idx], 0.25))/1.34

n = 25 # sample size
M=1000 # simulation size
B=1000 # bootstrap size

# Data structure 1
cover.mat <- matrix(NA, nrow = M, ncol = 3)
colnames(cover.mat) <- c('Normal', 'percentile', 'pivotal')

# data structure 2
len.mat <- matrix(NA,nrow = M, ncol = 3)
colnames(len.mat) <- c('Normal', 'percentile', 'pivotal')

# Main loop
set.seed(1231)
for(m in 1:M){
  # sample from a t-student with df=3
  x.sim <- rt(n,df=3)
  # plug-in estimator
  iqr.hat <- iqr(x.sim)
  # bootstrap
  t.boot <- boot(x.sim, iqr, B)
  
  # Calculate CI's
  normal <- c(iqr.hat - 2*sd(t.boot$t), iqr.hat + 2*sd(t.boot$t) )
  percen <- c(quantile(t.boot$t, 0.025), quantile(t.boot$t, 0.975))
  pivot <- c(2*iqr.hat - quantile(t.boot$t, 0.975),2*iqr.hat - quantile(t.boot$t, 0.025) )
  
  # Check if they cover or not
  cover.mat[m, 1] <- (iqr.true >= normal[1]) & (iqr.true <= normal[2])
  cover.mat[m, 2] <- (iqr.true >= percen[1]) & (iqr.true <= percen[2])
  cover.mat[m, 3] <- (iqr.true >= pivot[1]) & (iqr.true <= pivot[2])
  # store the length
  len.mat[m, 1] <- normal[2] - normal[1]
  len.mat[m, 2] <- percen[2] - percen[1]
  len.mat[m, 3] <- pivot[2] - pivot[1]
}

save(file = 'simulation_data.RData', cover.mat, len.mat)
#load('simulation_data.RData')
colMeans(cover.mat)

colMeans(len.mat)


len.df <- data.frame(Length = c(len.mat),
                     Type = factor(c(rep('Normal', M), rep('Percentile'), rep('Pivotal', M))))





# normal approxiamtion of the binomial distribution

n = 323
p = 0.35
M = 1e4

binom.sample <- rbinom(M, size = n, prob = p)

hist(scale(binom.sample), breaks = 25, probability = T, col= rgb(.5,0,.5,.2),
     main = 'Sample from a Bin(n,p')
curve(dnorm(x), add=T, col=6, lwd=3)



# Wald test for poisson

mu0 <- 1
n <- 20
a <- 0.05
M <- 1000
z.a <- qnorm(1 - a/2)

set.seed(123)

X <- matrix(rpois(n*M, lambda = mu0), nrow = M, ncol = n)

# MLE & SE

mle <- rowMeans(X)
se <- sqrt(mle/n)

# Wald test
W <- (mle - mu0)/se

# frequency of rejection
mean(abs(W)>z.a)




# Exercise 2

# compute the p-value of a t-test (comparison of variable)
MyStudent <- function(response, group){
  z = t.test(response ~ group, var.equal = T)
  z = z$p.value
  return(z)
}


require(mvtnorm, quietly = T)

n=30
m=500

set.seed(123)

dta0 =rmvnorm(2*n, mean = rep(0,m))

gr <- factor(rep(c(1,2), c(n,n)))

pval <- apply(dta0, 2, MyStudent, group = gr)

hist(pval, nclass = 10, col = 'darkgray')
