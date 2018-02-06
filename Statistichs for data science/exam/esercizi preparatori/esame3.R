# Question C

#Data
xx <- c(0.65, 0.48, 0.08, 0.11, 0.8, 0.67, 0.38, 0.63, 0.1, 0.03)
#MLE for psi
psi.ml <- function(x) - length(x)/sum(log(1- x))
round(psi.ml(xx),2)

#MLE for tau
tau.ml <- function(x) 1/(1 + psi.ml(x))
round(tau.ml(xx),2)

#Wald test

tau0 <- 1
#SE for tau 
se.tau <- function(x) (1/((1 + psi.ml(x))^2)) * ( abs(psi.ml(x))/sqrt(length(x)))
#Wald statistics
abs((tau.ml(xx) - tau0)/se.tau(xx))


# Likelihood function
L <- function(psi, x){
  n <- length(x)
  s <- prod((1 -x))
  out <- psi^n * s^(psi - 1)
  return(out)
}


# MLE
mle <- psi.ml(xx)

#Plot
curve(L(x, xx), from=0, to = 4, lwd = 4, col='orchid',
      xlab = expression(psi), ylab = 'Likelihood function')
segments(mle, 0, mle, L(mle, xx), lty = 3, col = 'pink2')
text(mle, 0.05, bquote(hat(psi)== .(round(mle, 2))), pos = 2, cex = .9)
grid()


# Parametric Boostrap
n <- length(xx)
B <- 1000
tau.boot  <- rep(NA, B)
set.seed(123)
for(b in 1:B){
  x.boot <- rbeta(n, 1, mle)
  tau.boot[b] <- 1/ (1 + psi.ml(x.boot))
}

# Bias via Parametric Bootstrap
cat('--- PBoot / Bias ---\n')
round(mean(tau.boot) - tau.ml(xx), 3)

# SE via Parametric Bootstrap
cat('--- PBoot / SE ---\n')
round(sd(tau.boot), 3)

# SE via Asymptotics
cat('--- Asymptotic SE ---\n')
round(se.tau(xx), 3)

# Parametric Bootstrap CI
cat('--- PBoot CI(0.9) --\n')
round(c(lower = tau.ml(xx) - qnorm(1 - 0.1/2)*sd(tau.boot),
        upper = tau.ml(xx) + qnorm(1 - 0.1/2)*sd(tau.boot)),3)


# Non-parametric Bootstrap
n <- length(xx)
B <- 1000
tau.nboot  <- rep(NA, B)
set.seed(123)
for(b in 1:B){
  idx <- sample(1:n, replace = T)
  x.nboot <- xx[idx]
  tau.nboot[b] <- 1/ (1 + psi.ml(x.nboot))
}

# Bias via Parametric Bootstrap
cat('--- NBoot / Bias ---\n')
round(mean(tau.nboot) - tau.ml(xx), 3)

# SE via Parametric Bootstrap
cat('--- NBoot / SE ---\n')
round(sd(tau.nboot), 3)

