# Stat4-DS2+CS

# March 21, 2017

# (last update: March 28, 2017)



# Code for A/R implementation in the particulare case:

# - the auxiliary is a uniform
# - the target is bounded on [0,1] (it is a Beta distribution)


curve(dunif(x)*3,0,1,xlab="x",ylab=expression(f[X](x)),ylim=c(0,4),lwd=2)
curve(dbeta(x,2,4),add=TRUE,col="red",lwd=2)

text(0.8,3.5,labels=expression(k~f[U](x)))
text(0.8,0.7,labels=expression(f[X](x)),col="red")

legend(x="topleft",lty=1,lwd=2.4,col=c("red","black"),legend=c("target density","bounding function"))
title(main="A/R")


### SIMULATION 

ef=function(x){
  dbeta(x,2,4)
}

k=3

n_sim_aux=10000

Y=rep(NA,n_sim_aux)
E=rep(NA,n_sim_aux)
for(i in 1:n_sim_aux){
  Y[i]=runif(1)
  E[i]=rbinom(1,size=1,prob=ef(Y[i])/k)
}

hist(Y,prob=TRUE)
hist(Y[E==1],prob=TRUE)



curve(dbeta(x,2,4),col="blue",lwd=2,add=TRUE)



### attempt to provide a general function

q=function(x){
  dunif(x)
}

draw_from_q=function(n){
  runif(n)
}

f=function(x){
  dbeta(x,2,4)
}


AR=function(dtarget,dauxiliary,rauxiliary,k){
  
  count=0
  E=0
  
  while(E==0){
    candidate = rauxiliary(1)
    acc_prob=dtarget(candidate)/(k*dauxiliary(candidate))
    E = sample(c(1,0),prob=c(acc_prob, 1-acc_prob),size=1)
    count=count+1
  }
  
  return(list(draw=candidate,computational_effort=count))
  
}


AR(dtarget=q,dauxiliary=q,rauxiliary=draw_from_q,k)

mcsize=1000
draw_vec=rep(NA,mcsize)
effort_vec=rep(NA,mcsize)

for(i in 1:mcsize){
  
  DD=AR(dtarget=f,dauxiliary=q,rauxiliary=draw_from_q,k=3)
  draw_vec[i] = DD$draw
  effort_vec[i] = DD$computational_effort
  
}

hist(draw_vec,freq=FALSE)
curve(f(x),add=TRUE)

plot(prop.table(table(effort_vec)),ylim=c(0,1),pch=16,col="red")
points(1:20,dgeom(0:19,prob=mean(E)))
