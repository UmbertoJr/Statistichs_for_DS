
#  question C -------------------------------------------------------------

## Plot the likelihood

#Data
xx <- c(-1,-1,0)
theta.hat <- function(x) mean(abs(x))
mle <- theta.hat(xx)

#Likelihood function
L <- function(theta, x){
  n <- length(x)
  jnk <- sum(abs(x))
  out <- ((theta/2)^jnk)*((1 - theta)^(n - jnk))
  return(out)
}

#Plot
curve(L(x,xx), from= 0, to=1, lwd=4, col='orchid', xlab = expression(theta), ylab = 'Likelihood function')
segments(mle, 0, mle, L(mle,xx), lty= 3, col='green')
text(mle,0.03, bquote(hat(theta)== .(round(mle,2))), pos=2, cex = .8)
grid()


## Write a function my.pmf(x, theta) that implements p(x|theta), and then use it together 
## with sample() to build a parametric bootstrap routine in order to evaluate bias and standard 
## error of MLE (based on MLE). Finally, find a 90% confidence interval for theta.

#Sampling from the disccrete model via <sample>
pmf <- function(x, theta){
  if((theta>1)||(theta<0)) stop("the parameter <theta> must be in (0,1)")
  support <- c(-1,0,1)
  out <- ((theta/2)^abs(x))*((1 - theta)^(1 - abs(x)))*(x %in% support)
  return(out)
}

# Parametric Boostrap
n <- length(xx)
B <- 1000
t.boot <- rep(NA, B)
set.seed(123)
for(b in 1:B){
  x.boot <- sample(c(-1,0,1), n, replace = T, prob = c(pmf(-1, mle),pmf(0, mle), pmf(1, mle)))
  t.boot[b] <- theta.hat(x.boot)
}

#Bias via Parametric Bootstrap
cat('---PBoot / Bias ----\n')
round(mean(t.boot) - mle, 3)

#SE via Parametric Bootstrap
cat('---PBoot / SE ----\n')
round(sd(t.boot), 3)

#Parametric bootstrap CI
cat('---PBoot CI(0.90) ----\n')
round(c(lower= mle - qnorm(1 - 0.1/2)*sd(t.boot),
        upper= mle + qnorm(1 - 0.1/2)*sd(t.boot)), 3)


#Sample
n <- 100
t.true <- 0.7 
dta <- sample(c(-1,0,1), n, replace = T, prob = c(pmf(-1, t.true),pmf(0, t.true), pmf(1, t.true)))
mle <- theta.hat(dta)


#Parametric bootstrap
B = 1000
t.boot <- rep(NA, B)
for(b in 1:B){
  x.boot <- sample(c(-1,0,1), n, replace = T, prob = c(pmf(-1, mle),pmf(0, mle), pmf(1, mle)))
  t.boot[b] <- theta.hat(x.boot)
}

#Bias via Parametric Bootstrap
cat('---PBoot / Bias ----\n')
round(mean(t.boot) - mle, 3)

#SE via Parametric Bootstrap
cat('---PBoot / SE ----\n')
round(sd(t.boot), 3)

#Parametric bootstrap CI
cat('---PBoot CI(0.90) ----\n')
round(c(lower= mle - qnorm(1 - 0.1/2)*sd(t.boot),
        upper= mle + qnorm(1 - 0.1/2)*sd(t.boot)), 3)


## Bonus

# Data
xx <- c(-1, -1, 0)
mle <- theta.hat(xx)
se.hat <- function(x) sqrt((theta.hat(x))*(1 - theta.hat(x))/length(x))


#Bias via Parametric Bootstrap
cat('---MLE ----\n')
round(mle, 3)

#SE via Parametric Bootstrap
cat('---SE ----\n')
round(se.hat(xx), 3)

#Parametric bootstrap CI
cat('---PBoot CI(0.90) ----\n')
round(c(lower= mle - qnorm(1 - 0.1/2)*se.hat(xx),
        upper= mle + qnorm(1 - 0.1/2)*se.hat(xx)), 3)

