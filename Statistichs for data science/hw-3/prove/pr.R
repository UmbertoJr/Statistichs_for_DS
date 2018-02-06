l = runif(100,1,3)
?runif
l <-seq(from = -1, to = 1,length.out = 1000)
plot(l,dunif(l,-1,1))
?dunif
l

set.seed(123)
M<-1000
z.samp = rnorm(M)
hist(z.samp,probability = T, breaks = 100)
curve(dnorm(x),col='red',add=T)


df.try=4
w.samp= rchisq(M, df = df.try)
hist(w.samp,probability = T,breaks = 100)
curve(dchisq(x, df=df.try),col= 'green',add=T)
lines(density(w.samp), col='blue')


t.samp = z.samp/sqrt(w.samp/df.try)

hist(t.samp, probability = T, breaks = 100)
lines(density(t.samp), col='blue')
curve(dt(x, df = df.try, ncp = 0),col='green', add=T)
curve(dnorm(x),col='red',add=T)

h=rpois(100,3)
hist(h, prob=T)
lines(density(h),col='blue')
plot(x,dpois(x,3),col='red',add = T)
?dpois
x=1:10
