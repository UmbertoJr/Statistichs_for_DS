# Question B Prob

#Check it's a density
integrate( function(x) 1/(2 *sqrt(1 - x)), 0,1)
integrate(function(y) 1/(2*y*sqrt(2 - log(y))), exp(1), exp(2))

#Take a look
par(mfrow=c(1,2))
curve(1/(2 *sqrt(1 - x)), 0,1, 
      col = 'purple', lwd= 4,
      main = 'density (before)', xlab = 'x', ylab = expression(f[X](x)))
curve(1/(2*x*sqrt(2 - log(x))), exp(1), exp(2), 
      col = 'lavender', lwd= 4,
      main = 'density (after)', xlab = 'y', ylab = expression(f[Y](y)))
