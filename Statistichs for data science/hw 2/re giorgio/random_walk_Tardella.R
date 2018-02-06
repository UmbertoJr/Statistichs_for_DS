# Last update 2017.04.10

set.seed(123)

# Very basic example of a Markov chain in discrete time: RANDOM WALK

# Let us consider the process X_t of 
# accumulation of wins or losses of one better

# each (discrete) time t=1,2,...

# you bet and 
# win    +1 with probability p 
# loose  -1 with probability 1-p 


########################################
### SIMULATION OF ONE SAMPLE PATH
########################################

p <- 0.7

index_set <-1:10

X <- seq_along(index_set)

X[1] <- 0 # INITIAL STATE (non random, although it need to to be)

for(t in index_set[-length(index_set)]){
  X[t+1] <- X[t] + sample(size=1,c(-1,1),prob=c(1-p,p))
}


########################################
### GRAPHICAL REPRESENTATION
########################################

# pdf(file="Plot-Sample-Path-random-walk.pdf",height = 5, width=9)
pdf(file="Plot-markov-present-past-future.pdf",height = 5, width=9)
plot(index_set,X,ylim=c(X[1]-max(index_set),X[1]+max(index_set)),type="b",pch=16,xlab="t",ylab=expression(X[t]))
title(main="Random walk - sample path")

highlight_present <- TRUE
if(highlight_present){
  eps=0.2
  rect(6-eps,-length(index_set)-2,6+eps,length(index_set)+2,col="yellow")
  points(index_set,X,ylim=c(X[1]-max(index_set),X[1]+max(index_set)),type="b",pch=16)
  
  text(6,-max(index_set),labels="present",srt=90,adj=0,cex=0.7)
  text(3,-max(index_set)+2,labels="past",srt=0,adj=0,cex=0.8)
  text(7.9,-max(index_set)+2,labels="future",srt=0,adj=0,cex=0.8)
  
}

newtemp <- c(0)

for(i in index_set){
  
  points(rep(i,length(newtemp)),newtemp,col="gray")
  rr=do.call("seq",as.list(range(newtemp)))
  rr=setdiff(rr,newtemp)
  points(rep(i,length(rr)),rr,col="gray",cex=0.3)
  newtemp=unique(c(sapply(newtemp,function(x) x + c(-1,+1))) )
  
}

points(1,X[1],pch=16,col="red")

dev.off()

