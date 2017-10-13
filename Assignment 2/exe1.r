#### Exe 1 ####
# E[(1-U^2)^(1/2)] with control variables
rm(list=ls())
# Z = U^2 vs Z
set.seed(1)
#do pilot simulation
p=100
Z = c(0,0)
cov_YZ = rep(0,2)
mean_Y=0
mean_Z=rep(0,2)
var_Z = rep(0,2)
for(i in 1:p){
  U = runif(1)
  Z = c(U^2,U)
  Y = (1-U^2)^(1/2)
  mean_Z = mean_Z + Z
  var_Z = var_Z + Z^2
  mean_Y = mean_Y + Y
  cov_YZ = cov_YZ + Y*Z
}
mean_Y = mean_Y/p
mean_Z = mean_Z/p
var_Z = (var_Z - p*mean_Z^2)/(p-1)
cov_YZ = (cov_YZ-p*mean_Z*mean_Y)/(p-1)

c = -cov_YZ/var_Z

# main simulation
n=1000
mean_V = 0
var_V = 0
for(i in 1:n){
  U = runif(1)
  Z = c(U^2,U)
  Y = (1-U^2)^(1/2)
  V = Y + c*(Z-mean_Z)
  mean_V = mean_V + V
  var_V = var_V + V^2
}

theta_hat = mean_V/n
var_V = (var_V - n*theta_hat^2)/(n-1)
a = 0.05
z_a = qnorm(1-a/2)
CI = matrix(c(theta_hat-z_a*sqrt(var_V/n),theta_hat+z_a*sqrt(var_V/n)),ncol=2)
