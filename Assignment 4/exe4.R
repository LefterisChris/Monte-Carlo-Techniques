#### exe 4 ####
rm(list=ls())

set.seed(1)
# standard estimator
N = 10000
X = runif(N)
Y = exp(X^2)
theta_0 = mean(Y)
var_0 = var(Y)

set.seed(1)
# stratified method
M = 10
Xj = rep(0,M)
tmp = rep(0,M)
for(i in 1:M){
  U = runif(N/M)/10 + 0.1*(i-1)
  Z = exp(U^2)
  Xj[i] = mean(Z)
  tmp[i] = var(Z)
}
theta_st = mean(Xj)
var_st = sum(tmp)/100/1000

opt_nj = round((tmp/10)/sum(tmp/10)*N)
Xj_opt = rep(0,M)
varj_opt = rep(0,M)
for(i in 1:M){
  U = runif(opt_nj[i])/10 + 0.1*(i-1)
  Z = exp(U^2)
  Xj_opt[i] = mean(Z)
  varj_opt[i] = var(Z)
}
theta_st_opt = mean(Xj_opt)
var_st_opt = sum(varj_opt/100/opt_nj)

set.seed(1)
theta_2 = 0
var_2 = 0
Nj = N/M
for(j in 1:M){
  sum = 0
  sum2 = 0
  for(i in 1:Nj){
    U = runif(1)/M+0.1*(j-1)
    Y = exp(U^2)
    sum = sum + Y
    sum2 = sum2 + Y^2
  }
  theta_j = sum/Nj
  var_j = (sum2-Nj*theta_j^2)/(Nj-1)
  theta_2 = theta_2+theta_j/10
  var_2 = var_2 + var_j/Nj/100
}