#### Exe 5 ####
library(Matrix)
library(matrixStats)
rm(list=ls())
Maturity =3
N = 10000
P = 1000
S0 = 100
r=0.04
sigma_low = 0.25
sigma_high = 0.75
lambda_low = 1
lambda_high = 3
K = c(50,100,150)

St = rep(S0,N)
Yt = rep(0,N)
Yt2 = rep(0,N)

for(i in 1:N){
  t = 0
  curr_sigma = sigma_low
  curr_lambda = lambda_low
  new_sigma = sigma_high
  new_lambda = lambda_high
  while (t < Maturity){
    U = runif(1)
    t_jump = -log(U)/curr_lambda
    if (t+t_jump <= Maturity){
      dlt = t_jump
      t = t + t_jump
    }
    else{
      dlt = Maturity - t
      t = Maturity
    }
    Z = rnorm(1)
    Yt[i] = Yt[i]+((r-sigma_low^2/2)*dlt+sigma_low*sqrt(dlt)*Z)
    Yt2[i] = Yt2[i]+((r-sigma_low^2/2)*dlt+sigma_low*sqrt(dlt)*Z)^2
    
    St[i] = St[i]*exp((r-curr_sigma^2/2)*dlt+curr_sigma*sqrt(dlt)*Z)
    
    tmp = curr_sigma
    curr_sigma = new_sigma
    new_sigma = tmp
    tmp = curr_lambda
    curr_lambda = new_lambda
    new_lambda = tmp
  }
}
Yt_exp = exp(Yt)
CTRL = matrix(0,nrow = N,ncol = length(K))
for(i in 1:length(K)){
  CTRL[,i] = rowMaxs(matrix(c(rep(0,N),S0*Yt_exp-K[i]),nrow = N))
}

mu_St = mean(St)
mu_Y = mean(Yt)
mu_Y2 = mean(Yt2)
mu_Yexp = mean(Yt_exp)
mu_CTRL = colMeans(CTRL)

# Initial estimator
C = matrix(0,nrow = N,ncol = length(K))
for(i in 1:length(K)){
  tmp = matrix(c(rep(0,N),St-K[i]),ncol=2)
  C[,i] = exp(-r*Maturity)*rowMaxs(tmp)
}

c1 = as.vector(-cov(St,C)/var(St))
c2 = as.vector(-cov(Yt,C)/var(Yt))
c3 = as.vector(-cov(Yt2,C)/var(Yt2))
c4 = matrix(0,nrow = 2,ncol = length(K))
for(i in 1:length(K)){
  lg = lm(C[,i]~Yt+Yt2)
  c4[,i] = -lg$coefficients[2:3]
}
c5 = as.vector(-cov(Yt_exp,C)/var(Yt_exp))
c6 =matrix(0,nrow=2,ncol = length(K))
for(i in 1:length(K)){
  lg2 = lm(C[,i]~St+Yt)
  c6[,i]=-lg2$coefficients[2:3]
}
c7 = -diag(cov(CTRL,C))/diag(var(CTRL))
c8 = matrix(0,nrow = 2,ncol = length(K))
for(i in 1:length(K)){
  lg = lm(C[,i]~St+CTRL[,i])
  c8[,i] = -lg$coefficients[2:3]
}
c9 = matrix(0,nrow = 3, ncol = length(K))
for(i in 1:length(K)){
  lg = lm(C[,i]~St+Yt+CTRL[,i])
  c9[,i] = -lg$coefficients[2:4]
}

C = matrix(0,nrow = N,ncol = length(K))
C1 = matrix(0,nrow = N,ncol = length(K))
C2 = matrix(0,nrow = N,ncol = length(K))
C3 = matrix(0,nrow = N,ncol = length(K))
C4 = matrix(0,nrow = N,ncol = length(K))
C5 = matrix(0,nrow = N,ncol = length(K))
C6 = matrix(0,nrow = N,ncol = length(K))
C7 = matrix(0,nrow = N,ncol = length(K))
C8 = matrix(0,nrow = N,ncol = length(K))
C9 = matrix(0,nrow = N,ncol = length(K))
for(i in 1:length(K)){
  tmp = matrix(c(rep(0,N),St-K[i]),ncol=2)
  C[,i] = exp(-r*Maturity)*rowMaxs(tmp)
  C1[,i] = C[,i]+c1[i]*(St-mu_St)
  C2[,i] = C[,i]+c2[i]*(Yt-mu_Y)
  C3[,i] = C[,i]+c3[i]*(Yt2-mu_Y2)
  C4[,i] = C[,i] + c4[1,i]*(Yt-mu_Y) + c4[2,i]*(Yt2-mu_Y2)
  C5[,i] = C[,i] + c5[i]*(Yt_exp - mu_Yexp)
  C6[,i] = C[,i] + c6[1,i]*(St-mu_St) + c6[2,i]*(Yt-mu_Y)
  C7[,i] = C[,i] + c7[i]*(CTRL[,i]-mu_CTRL[i])
  C8[,i] = C[,i] + c8[1,i]*(St-mu_St) + c8[2,i]*(CTRL[,i]-mu_CTRL[i])
  C9[,i] = C[,i] + c9[1,i]*(St-mu_St) + c9[2,i]*(Yt-mu_Y) + c9[3,i]*(CTRL[,i]-mu_CTRL[i])
}

cat("Prices:",colMeans(C),"\n")
cat("Regular variances:",colVars(C),"\n")
cat("Prices with control variable St:",colMeans(C1),"\n")
cat("Control variable St:",colVars(C1),"\n")
cat("Prices with control variable Yt:",colMeans(C2),"\n")
cat("Control variable Yt:",colVars(C2),"\n")
cat("Prices with control variable Yt^2:",colMeans(C3),"\n")
cat("Control variable Yt^2:",colVars(C3),"\n")
cat("Prices with control variable Yt+Yt^2:",colMeans(C4),"\n")
cat("Control variable Yt+Yt^2:",colVars(C4),"\n")
cat("Prices with control variable exp(Yt):",colMeans(C5),"\n")
cat("Control variable exp(Yt):",colVars(C5),"\n")
cat("Prices with control variable St+Yt:",colMeans(C6),"\n")
cat("Control variable St+Yt:",colVars(C6),"\n")
cat("Prices with control variable CTRL:",colMeans(C7),"\n")
cat("Control variable CTRL:",colVars(C7),"\n")
cat("Prices with control variable St+CTRL:",colMeans(C8),"\n")
cat("Control variable St+CTRL:",colVars(C8),"\n")
cat("Prices with control variable St+CTRL:",colMeans(C9),"\n")
cat("Control variable St+CTRL:",colVars(C9),"\n")
