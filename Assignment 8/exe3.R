#### exe 3 ####
# this function has argument a matrix MxN where M/2 where the Markov chains 
# and each had 2*N samples from the distribution. So each MC produces two MC's
# with N samples each
Gelman_Rubin_diagnostic <- function(param){
  dims = dim(param)
  chain_means = rowMeans(param)
  grand_mean = mean(chain_means)
  # between-sequence variance
  B = dims[2]/(dims[1]-1)*sum((chain_means-grand_mean)^2)
  # within-sequence variance
  tmp = 0.0
  for(i in 1:dims[1]){
    tmp = tmp + sum((param[i,]-chain_means[i])^2)/(dims[2]-1)
  }
  W = tmp/dims[1]
  # var hat
  var_hat = (dims[2]-1)/dims[2]*W+B/dims[2]
  # R statistical summary of MCMC that diagnoses convergenve. 
  # Good values close to 1
  R = sqrt(var_hat/W)
  return(R)
}

p_party = 0.2
p_not_motivated = 0.4
p_angry = c(0.95,0.5)
p_headache = c(0.2,0.9)
p_under_perf = c(0.999,0.9,0.9,0.01)

calc_joint_dist <- function(theta){
  # find the pmf for under performing
  p_U=(p_under_perf^theta[3]*(1-p_under_perf)^(1-theta[3]))
  # calculate the joint distribution
  joint_dist = ((p_party^theta[1] * (1 - p_party)^(1 - theta[1])) * 
        (p_not_motivated^theta[2]*(1 - p_not_motivated)^(1 - theta[2])) *
        (p_headache[1]^(1 - theta[1])*p_headache[2]^theta[1]) * 
        (p_U[1]^(theta[1]*theta[2]))*(p_U[2]^(theta[1]*(1-theta[2]))) * 
        (p_U[3]^((1-theta[1])*theta[2]))*(p_U[4]^((1-theta[1])*(1-theta[2]))) * 
        (p_angry[1]^theta[3] * p_angry[2]^(1-theta[3])))
  return(joint_dist)
}

Metropolis_Hastings <- function(theta,N0,N) {
  # returns theta parameter vector
  Tot_iter = N0 + N
  vec_theta = matrix(0,nrow=N,ncol=3)
  for(iter in 1:Tot_iter){
    new_theta = theta
    # randomly pick one of {P,D,U}
    U = runif(1)
    if(U<1/3){
      # we choose P - flip it
      new_theta[1] = !theta[1]
    } else if(U<2/3){
      # we choose D - flip it
      new_theta[2] = !theta[2]
    } else{
      # we choose U - flip it
      new_theta[3] = !theta[3]
    }
    # calculate the target dist with the old generated parameter vector
    p_old = calc_joint_dist(theta)
    # calculate the target distribution with the new generated parameter vector
    p_new = calc_joint_dist(new_theta)
    # calculate the acceptance probability - the jumping distribution equals 1/3 so it cancels out
    r = min(p_new/p_old,1)
    # test acceptance - rejection
    U = runif(1)
    if(U<r){
      theta = new_theta
    }
    if(iter > N0){
      # remove the burn in samples
      vec_theta[iter-N0,] = theta  
    }
  }
  return(vec_theta)
}

# initialize parameter vector theta = (P,D,U)
theta = c(0,0,0)

N = 50000
N0 = N

MH_theta = Metropolis_Hastings(theta,N0,N)
means=colMeans(MH_theta)
print(means[1])

init_thetas = matrix(c(0,0,0,0,1,0,1,0,1,1,1,1),nrow=4,byrow = T)
samples = matrix(0,nrow = 8,ncol=N/2)
for(i in 1:4){
  # theta = sample(0:1,3,replace = T)
  theta = init_thetas[i,]
  MH_theta = Metropolis_Hastings(theta,N0,N)
  samples[i,] = MH_theta[1:(N/2),1]
  samples[i+4,] = MH_theta[(N/2+1):N,1]
}
diag=Gelman_Rubin_diagnostic(samples)
print(diag)

# part b - Gibbs sampling

Gibbs_sampling <- function(theta,N0,N){
  Tot_iter = N0 + N
  K = length(theta)
  samples = matrix(0,nrow = N,ncol=K)
  for(iter in 1:Tot_iter){
    for(var in 1:K){
      # sample from the conditional distribution for every parameter
      tmp1 = tmp2 = theta
      tmp1[var] = 0
      tmp2[var] = 1
      # calculate the conditional distribution with the parameter equal to 0
      p0 = calc_joint_dist(tmp1)/(calc_joint_dist(tmp1)+calc_joint_dist(tmp2))
      U = runif(1)
      if (U<p0){
        theta[var] = 0
      }
      else{
        theta[var] = 1
      }
    }
    if(iter > N0){
      # burn-in samples are ommited
      samples[iter-N0,] = theta
    }
  }
  return(samples)
}

# initialize the parameter vector - {P,D,U}
theta = c(0,0,0)
N = 50000
N0 = N
G_theta = Gibbs_sampling(theta,N,N0)
colMeans(G_theta)

