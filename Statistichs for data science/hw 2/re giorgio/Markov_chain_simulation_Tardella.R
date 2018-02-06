#########################################################
# Markov Chain simulation on a discrete state space
# with only two states
# S={1,2}
#########################################################

S=c(1,2)        # discrete state space


alpha<-0.2      # set up the TPM (Transition Probability Matrix)
beta<-0.5
tpm<-matrix(c(1-alpha,alpha,beta,1-beta),nrow=2,byrow=T)

tpm


x1<-2 ## starting value of the chain
## [from a degenerate distribution with mass 1 at state s=2]

nsample<-10
chain<-rep(NA,nsample+1) # vector that will hold
# all the simulated values
chain[1]<-x1             # starting value x1 assigned to chain[1]
for(t in 1:nsample){
  chain[t+1]<-sample(S,size=1,prob=tpm[chain[t],])
}

S
plot(chain,ylim=c(0,4))
plot(chain,ylim=c(0,4),type="b")
table(chain)

prop.table(table(chain))

mean(chain)


### Let us start providing some evidence of stability of the behaviour of the initial stream of one sample path of the Markov chain of suitable length

nsample<-10000 #  length of the initial part of the sample path 
#  of our Markov chain (stochastic process) 

chain<-rep(NA,nsample+1) # vector that will hold
# all the simulated values

# Let us start from x1=1

x1=1
chain[1]<-x1             # starting value x1 assigned to chain[1]
for(t in 1:nsample){
  chain[t+1]<-sample(S,size=1,prob=tpm[chain[t],])
}

mean(chain)       # our first example of empirical mean [h(x)=x]
mean(chain==2)    # our second example of empirical mean [h(x)=I_{1}(x)]
table(chain)      
prop.table(table(chain))



# Let us start the chain from another value x1=2

x1=2
chain[1]<-x1             # starting value x1 assigned to chain[1]
for(t in 1:nsample){
  chain[t+1]<-sample(S,size=1,prob=tpm[chain[t],])
}

mean(chain)       # our first example of empirical mean [h(x)=x]
mean(chain==2)    # our second example of empirical mean [h(x)=I_{1}(x)]
table(chain)      
prop.table(table(chain))

# we experienced **approximately** the same empirical mean no matter what is the initial POINT we start the chain from!! [initial distribution can be regarded as and is in fact degenerate at that point!!]

# Let us look at the stabilization of the empirical mean as long as the number of sampled states increases

runningmeans=cumsum(chain)/(1:length(chain))
plot(1:length(chain),runningmeans)
title(main="stabilization of the running means")
runningmeans[length(runningmeans)]


# Let us look at the stabilization of the empirical mean of the indicator function of state s=1 (relative frequency) as long as the number of sampled states increases

state=1
runningmeans=cumsum(chain==state)/(1:length(chain))
plot(1:length(chain),runningmeans)
title(main="stabilization of the running relative frequencies \n of occurrence of state 1")
runningmeans[length(runningmeans)]

# now let us compare this value with the following entries 
# and try to interpret what is going on 

tpm
tpm%*%tpm
tpm%*%tpm%*%tpm
tpm%*%tpm%*%tpm%*%tpm
tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm
# ...
tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm
tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm%*%tpm








#########################################################
# Markov Chain simulation on a continuous state space
# AR(1) process [random walk]
#########################################################

# X_{t+1} = rho * X_{t+1} + B_{t+1}
# B_1,...,B_{t+1} i.i.d. N(0,1)

# ==> X_{t+1} ~ N(rho*X_{t},sd=1)

set.seed(123)
mcsize=10000
# Setup of the Markov transition kernel correponding 
# to the distribution of X_{t+1}|X_{t}=present

rho=1.002
rho=-0.2
rho=0.9

starting_point=present=3.9
chain=c(starting_point)

# initialization of the chain at time t=0

# loop for sequential updating from present to the next future time 
# transition kernel (updating) 
# made of a gaussian random draw where the gaussian distribution
# has a mean depending on the present state

print(present)

for(i in 1:mcsize){
  
  # future=rnorm(1,mean=0.2*present,sd=sqrt(1))
  # future=rnorm(1,mean=present)
  # future=rnorm(1,mean=1.05*present)
  # future=rnorm(1,mean=2*present)
  
  # alternative kernels depending on rho
  # rho=0.9
  
  future=rnorm(1,mean=rho*present,sd=sqrt(1))
  present=future
  chain=c(chain,present)
  # print(present)
  
}

# final print of the realization of the path up to time t (t=mcsize)
# (realization of the initial portion of the chain)

plot(chain)

plot(chain,type="b")

iterationvec=seq(1,length(chain))

plot(iterationvec,chain,type="l",main="Trace plot", sub="AR(1) - Random walk",xlab="t",ylab="chain state at time t")

hist(chain,freq=F)
curve(dnorm(x,sd=1/sqrt(1-rho^2)),col="red",add=T)

# 

acf(chain)
acf(chain)$acf[2]

# This continuous case example teaches us that things can go well eventually when n increases but not necessarily .... depending on ... what?

