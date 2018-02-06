
# Question C

#Ex-1

dmy <- function(x, psi) (x/psi^2)*exp( -x^2/(2*psi^2))*(x>0)

#Plot

psi.seq <- c(1,5,10)
curve(dmy(x, psi.seq[1]), from = 0, to = 30, n = 300, 
      lwd = 4, col = 'purple', ylab = expression(f[psi](x)))
curve(dmy(x, psi.seq[2]),lwd = 4, col = 'blue', add = T)
curve(dmy(x, psi.seq[3]),lwd = 4, col = 'pink', add = T)


#Ex-2

rmy <- function(n , psi) psi* sqrt(-2* log(1 - runif(n)))

#Try it out
set.seed(123)
M = 1000
par(mfrow= c(1,3))

# psi = 1
hist(rmy(M, psi.seq[1]), probability = T, border = 'white', col = 'pink',
     xlab = 'x', ylab = expression(f[psi](x)), main = bquote(psi == .(psi.seq[1])))
curve(dmy(x, psi.seq[1]), col= 'purple', add=T)

# psi = 5
hist(rmy(M, psi.seq[2]), probability = T, border = 'white', col = 'pink',
     xlab = 'x', ylab = expression(f[psi](x)), main = bquote(psi == .(psi.seq[2])))
curve(dmy(x, psi.seq[2]), col= 'purple', add=T)

# psi = 10
hist(rmy(M, psi.seq[3]), probability = T, border = 'white', col = 'pink',
     xlab = 'x', ylab = expression(f[psi](x)), main = bquote(psi == .(psi.seq[3])))
curve(dmy(x, psi.seq[3]), col= 'purple', add=T)

dev.off()

#Ex-3

#Setup
tau0 <- 15
psi0 <- sqrt(tau0 / (2 - (pi/2)))
n <- 100
a <- 0.05
M <- 1e3
z.a <- qnorm(1- a/2)

#Vectorize the simulation under the null hypothesis
set.seed(432)
X <- matrix(rmy(n*M, psi = psi0), nrow = M, ncol = n)

#MLE's & SE: functions
tau <- function(psi) (2 - pi/2)*(psi^2)
mle <- function(x) sqrt(sum(x^2)/(2*length(x)))
se.tau <- function(x) (abs(2*(2 -pi/2)*mle(x)))*(abs(mle(x))/(2*sqrt(length(x))))

#MLE's & SE: result
psi.mle.vec <- apply(X, 1, mle)
tau.mle.vec <- tau(psi.mle.vec)
tau.se.vec <- apply(X,1, se.tau)

# Test statistics
W <- abs((tau.mle.vec - tau0)/tau.se.vec)

#The frequency with which we *reject* the (true!) null should be around 0.05
mean(abs(W) > z.a)



