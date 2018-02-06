library(markovchain)

states <- c(1,2,3)
ByRow <- TRUE
transition_matrix <- matrix(data= c(0,1/2,1/2,
                                    5/8,1/8,1/4,
                                    2/3,1/3,0), byrow = ByRow,nrow=3)
Mc_states <- new('markovchain', states = states, byrow = ByRow, transitionMatrix = transition_matrix, name='Position')

initial_state <- c(1,0,0)



pi <- Mc_states^1000



plot(Mc_states)

n=100
after_n <- initial_state*(Mc_states^100)

after_n

  
states <- c(1,2,3)
transition_matrix <- matrix(data= c(0,1/2,1/2,
                                    5/8,1/8,1/4,
                                    2/3,1/3,0), byrow = TRUE ,nrow=3)




x1 <- 1

nsample<-1000
chain<-rep(NA,nsample+1) # vector that will hold
# all the simulated values
chain[1]<-x1             # starting value x1 assigned to chain[1]
for(t in 1:nsample){
  chain[t+1]<-sample(states,size=1,prob=transition_matrix[chain[t],])
}





t0 <- 1

ntimes <- 500
nsample<-1000
t1000 <-rep(NA,ntimes)

for(i in 1:ntimes){
  chain<-rep(NA,nsample+1) # vector that will hold
  chain[1]<-t0             # starting value x1 assigned to chain[1]
  for(t in 1:nsample){
    chain[t+1]<-sample(states,size=1,prob=transition_matrix[chain[t],])
  }
  t1000[i]<- chain[nsample+1]
}
