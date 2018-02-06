
# Evaluate the PMF of a Poisson with a very large parameter: < lambda = 1472 >
xx = 1300:1650
PMF.pois = dpois(xx , lambda = 1472)

# Plot  the Poisson(1472)
plot(xx, PMF.pois, type='h', lwd=.1,
     xlab = expression(lambda==1472), ylab = '',
     main = 'Poisson vs Normal')

#Add points
points(xx, PMF.pois, pch= '*', cex=1.5, col='red')

#Add the density function of a Normal
curve(dnorm(x, mean = 1472, sd = sqrt(1472)), add=T, lwd= 6, col=rgb(0,1,1,0.2))


### 27.2

pour.temp <- c(2543,2541,2544,2620,2560,2559,2562,2553,2552,2553)
mean(pour.temp)

#take a look to help 
t.test(pour.temp, mu = 2550, alternative = 'two.sided', conf.level = 0.99)


# Exercise 7
#data
x1 <- c(.225, .262, .217, .240, .230, .229, .235, .217)
x2 <- c(.209, .205, .196, .210, .202, .207, .224, .223, .220, .201)
y <- c(rep('Twain', length(x1)), rep('Snodgrass', length(x2)))
dta <- data.frame(Freq = c(x1, x2), Autor= y)
str(dta)

# Take a look
with(dta,
     boxplot(Freq ~ Autor, horizontal = T,
             main = 'frequency', col= c('purple', 'gold')))

aggregate(Freq ~ Autor, data= dta, FUN= var)

t.test(Freq ~ Autor,data = dta, var.equal=T)
