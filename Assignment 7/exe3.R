#### Exe 3 ####

Gelman_Rubin_diagnostic <- function(param) {
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


data = matrix(c(62,60,63,59,0,0,0,0,63,67,71,64,65,66,0,0,68,66,71,67,68,68,0,0,52,62,60,61,63,64,63,59),
              nrow=4,byrow = T)

J = dim(data)[1]
n = length(which(data!=0))

# number of chains
M = 8
# number of iterations
N = 1000
# number of burn-in iterations
N0 = 2*N
# total iterations
tot_iter = 2*N + N0

# define matrix for parameters with samples per chain
theta1_samples = matrix(0.0,nrow = M,ncol = N)
theta2_samples = matrix(0.0,nrow = M,ncol = N)
theta3_samples = matrix(0.0,nrow = M,ncol = N)
theta4_samples = matrix(0.0,nrow = M,ncol = N)
mu_samples = matrix(0.0,nrow = M,ncol = N)
sigma_samples = matrix(0.0,nrow = M,ncol = N)
tau_samples = matrix(0.0,nrow = M,ncol = N)

#initialize thetas with over-dispersed values
thetas = matrix(0,nrow=M/2,ncol=J)
for(i in 1:J){
  thetas[,i] = sample(data[i,which(data[i,]!=0)],size = 4)
}

for(chain in 1:(M/2)){
  # for each chain
  mu = mean(thetas[chain,])
  # parameters: theta1,theta2,theta3,theta4,mu,sigma,tau
  parameters = c(thetas[chain,],mu,0,0)
  
  # begin Gibbs sampling
  for(iter in 1:tot_iter){
    # for each iteration 2*N+N0
    # generate from conditional posterior distribution of sigma
    tmp = 0;
    for(i in 1:J){
      tmp = tmp + sum((data[i,which(data[i,]!=0)]-parameters[i])^2)/n  
    }
    sigma_2 = n*tmp/rchisq(1,df=n)
    # parameters[J+2] = sqrt(sigma_2)
    parameters[J+2] = sigma_2
    
    # generate from conditional posterior distribution of tau
    tmp = sum((parameters[1:J]-parameters[J+1])^2)/(J-1)
    tau_2 = (J-1)*tmp/rchisq(1,df=(J-1))
    # parameters[J+3] = sqrt(tau_2)
    parameters[J+3] = tau_2
    
    # generate from conditional posterior distribution of mu
    tmp = sum(parameters[1:J])/J
    # mu = rnorm(n=1,mean=tmp,sd = sqrt(parameters[J+3]^2/J))
    mu = rnorm(n=1,mean=tmp,sd = sqrt(parameters[J+3]/J))
    parameters[J+1] = mu
    
    for(k in 1:J){
      #generate from conditional posterior distribution of theta j
      idx = which(data[k,]!=0)
      # theta_j_hat = (mu/parameters[J+3]^2+length(idx)*mean(data[k,idx])/parameters[J+2]^2)/(1/parameters[J+3]^2+length(idx)/parameters[J+2]^2)
      # var_j_hat = 1/(1/parameters[J+3]^2+length(idx)/parameters[J+2]^2)
      theta_j_hat = (mu/parameters[J+3]+length(idx)*mean(data[k,idx])/parameters[J+2])/(1/parameters[J+3]+length(idx)/parameters[J+2])
      var_j_hat = 1/(1/parameters[J+3]+length(idx)/parameters[J+2])
      parameters[k] = rnorm(n=1,mean = theta_j_hat,sd = sqrt(var_j_hat))
    }
    if(iter > N0){
      # we discard the first N0 samples as burn-in samples
      if(iter <= (N0+N)){
        # split the chain into two each containing N samples.
        theta1_samples[chain,iter-N0] = parameters[1]
        theta2_samples[chain,iter-N0] = parameters[2]
        theta3_samples[chain,iter-N0] = parameters[3]
        theta4_samples[chain,iter-N0] = parameters[4]
        mu_samples[chain,iter-N0] = parameters[5]
        sigma_samples[chain,iter-N0] = parameters[6]
        tau_samples[chain,iter-N0] = parameters[7]
      }
      else{
        theta1_samples[chain+M/2,iter-N0-N] = parameters[1]
        theta2_samples[chain+M/2,iter-N0-N] = parameters[2]
        theta3_samples[chain+M/2,iter-N0-N] = parameters[3]
        theta4_samples[chain+M/2,iter-N0-N] = parameters[4]
        mu_samples[chain+M/2,iter-N0-N] = parameters[5]
        sigma_samples[chain+M/2,iter-N0-N] = parameters[6]
        tau_samples[chain+M/2,iter-N0-N] = parameters[7]
      }
    }
  }
}

theta1 = round(sort(colMeans(theta1_samples)),2)
theta2 = round(sort(colMeans(theta2_samples)),2)
theta3 = round(sort(colMeans(theta3_samples)),2)
theta4 = round(sort(colMeans(theta4_samples)),2)
mu = round(sort(colMeans(mu_samples)),2)
sigma = round(sort(sqrt(colMeans(sigma_samples))),2)
tau = round(sort(sqrt(colMeans(tau_samples))),2)
# run the diagnostic test
cat('Estimand\t\t Posterior Quantiles\t R\n')
cat('\t\t 2.5%  25%   50%   75%   97.5%\n')
cat('-------------------------------------------\n')
cat('theta1\t\t',theta1[floor(0.025*N)],theta1[floor(0.25*N)],theta1[floor(0.5*N)],theta1[floor(0.75*N)],theta1[floor(0.975*N)],round(Gelman_Rubin_diagnostic(theta1_samples),3),'\n')
cat('theta2\t\t',theta2[floor(0.025*N)],theta2[floor(0.25*N)],theta2[floor(0.5*N)],theta2[floor(0.75*N)],theta2[floor(0.975*N)],round(Gelman_Rubin_diagnostic(theta2_samples),3),'\n')
cat('theta3\t\t',theta3[floor(0.025*N)],theta3[floor(0.25*N)],theta3[floor(0.5*N)],theta3[floor(0.75*N)],theta3[floor(0.975*N)],round(Gelman_Rubin_diagnostic(theta3_samples),3),'\n')
cat('theta4\t\t',theta4[floor(0.025*N)],theta4[floor(0.25*N)],theta4[floor(0.5*N)],theta4[floor(0.75*N)],theta4[floor(0.975*N)],round(Gelman_Rubin_diagnostic(theta4_samples),3),'\n')
cat('mu    \t\t',mu[floor(0.025*N)],mu[floor(0.25*N)],mu[floor(0.5*N)],mu[floor(0.75*N)],mu[floor(0.975*N)],round(Gelman_Rubin_diagnostic(mu_samples),3),'\n')
cat('sigma \t\t',sigma[floor(0.025*N)],sigma[floor(0.25*N)],sigma[floor(0.5*N)],sigma[floor(0.75*N)],sigma[floor(0.975*N)],round(Gelman_Rubin_diagnostic(sigma_samples),3),'\n')
cat('tau   \t\t',tau[floor(0.025*N)],tau[floor(0.25*N)],tau[floor(0.5*N)],tau[floor(0.75*N)],tau[floor(0.975*N)],round(Gelman_Rubin_diagnostic(tau_samples),3),'\n')

