#### exe 1 ####
rm(list=ls())
library(MCMCpack)

n_A = 186
n_B = 38
n_AB = 13
n_O = 284

# initial values
probs = c(0.3,0.3,0.4)

N = 10000
N0 = N
Tot_iter = N + N0
samples = matrix(0,nrow=N,ncol=3)
for(iter in 1:Tot_iter){
  Z_A = rbinom(1,n_A,probs[1]^2/(probs[1]^2+2*probs[1]*probs[3]))
  Z_B = rbinom(1,n_B,probs[2]^2/(probs[2]^2+2*probs[2]*probs[3]))
  new_probs = rdirichlet(1,alpha = c(Z_A+n_A+n_AB+1,Z_B+n_B+n_AB+1,2*n_O+n_A+n_B-Z_A-Z_B+1))
  probs = new_probs
  if(iter > N0){
    samples[iter-N0,] = probs
  }
}
par(mfrow=c(1,3))
hist(samples[,1],freq = F,main = 'Histogram for p_A',xlab = 'p_A')
parhist(samples[,2],freq=F,main='Histogram for p_B',xlab = 'p_B')
hist(samples[,3],freq=F,main='Histogram for p_O',xlab = 'p_O')
