#### Exe 1 ####
# E[(1-U^2)^(1/2)] with control variables
rm(list=ls())
# Z = U^2 vs Z
set.seed(1)
#do pilot simulation
p=1000
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
n=10000
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

# method with linear regression 
# rm(list=ls())
# set.seed(1)
# tmp = runif(1000)
# y = (1-tmp^2)^(1/2)
# x1 = tmp^2
# x2 = tmp
# lg_YZ1=lm(y~x1)
# lg_YZ2=lm(y~x2)


#### Exe 2 ####
# with two covariates U, U^2
rm(list=ls())

# method with linear regression
set.seed(1)
U = runif(1000)
Y = (1-U^2)^(1/2)
Z = matrix(c(U,U^2),ncol = 2)
lg1 = lm(Y~Z)
c = -lg1$coefficients[-1]

n=10000
U = runif(n)
Y = (1-U^2)^(1/2)
Z = matrix(c(U,U^2),ncol=2)
mean_Z = c(mean(Z[,1]),mean(Z[,2]))
var_Z = c(var(Z[,1]),var(Z[,2]))
theta_hat = mean(Y+c[1]*(Z[,1]-mean_Z[1])+c[2]*(Z[,2]-mean_Z[2]))
var(Y+c[1]*(Z[,1]-mean_Z[1])+c[2]*(Z[,2]-mean_Z[2]))

#### Exe 3 ####
rm(list=ls())
set.seed(1)
P <-10000
n <-1000
m <-11
S0 <-100
r <- .05
T<-1
sigma <-.25 
dt = T/n
GBM<-matrix(0,nrow = P,ncol = n)
GBM[,1] <- rep(S0,P)
for(i in 1:(n-1)){
  GBM[,i+1] <- GBM[,i] * exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*rnorm(P)) 
}
K = c(90,100,110,120)

# group<-which(1:n %% floor(n/m)==0)
# group <- sample(1:n,m)
group <- round(seq(from=floor(n/m),to=n,length.out = m))

a = 0.99
z_a = qnorm((1+a)/2)

for(i in 1:length(K)){
  payoff <-0
  payoff_2 <-0
  for(j in 1:P){
    priceT <- sum(GBM[j,group])/m
    tmp = exp(-r*T)*max(0,priceT-K[i])
    payoff <- payoff + tmp
    payoff_2 <- payoff_2 + tmp^2 
  }
  mean_P <- payoff/P
  var_P <- (payoff_2 - P*mean_P^2)/(P-1)
  CI <- c(mean_P-z_a*sqrt(var_P/P),mean_P+z_a*sqrt(var_P/P))
  cat("For strike K=",K[i],":",mean_P,"with",100*a,"% confidence interval: [",CI[1],",",CI[2],"]\n")
}

# with control variable Z = SUM(S_{iT/m})/m
set.seed(1)
pilot=1000
GBM_pilot<-matrix(0,nrow = pilot,ncol = n)
GBM_pilot[,1] <- rep(S0,pilot)
for(i in 1:(n-1)){
  GBM_pilot[,i+1] <- GBM_pilot[,i] * exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*rnorm(pilot)) 
}
c = rep(0,length(K))
mean_Z = rep(0,length(K))
var_Z = rep(0,length(K))
cov_YZ=rep(0,length(K))
for(i in 1:length(K)){
  payoff = 0
  for(j in 1:pilot){
    priceT <- sum(GBM_pilot[j,group])/m
    tmp = exp(-r*T)*max(0,priceT-K[i])
    payoff <- payoff + tmp
    cov_YZ[i] = cov_YZ[i] + priceT*tmp 
    mean_Z[i] = mean_Z[i] + priceT
    var_Z[i] = var_Z[i] + priceT^2
  }
  mean_P <- payoff/pilot
  mean_Z[i] <- mean_Z[i]/pilot
  var_Z[i] = (var_Z[i]-pilot*mean_Z[i]^2)/(pilot-1)
  cov_YZ[i] <- (cov_YZ[i]-pilot*mean_Z[i]*mean_P)/(pilot-1)
  c[i] = -cov_YZ[i]/var_Z[i]
}

# GBM<-matrix(0,nrow = P,ncol = n)
# GBM[,1] <- rep(S0,P)
# for(i in 1:(n-1)){
#   GBM[,i+1] <- GBM[,i] * exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*rnorm(P))
# }

for(i in 1:length(K)){
  payoff <-0
  payoff_2 <-0
  for(j in 1:P){
    priceT <- sum(GBM[j,group])/m
    tmp = exp(-r*T)*max(0,priceT-K[i])
    V = tmp +c[i]*(priceT-mean_Z[i])
    payoff <- payoff + V
    payoff_2 <- payoff_2 + V^2 
  }
  mean_P <- payoff/P
  var_P <- (payoff_2 - P*mean_P^2)/(P-1)
  CI <- c(mean_P-z_a*sqrt(var_P/P),mean_P+z_a*sqrt(var_P/P))
  cat("For strike K=",K[i],":",mean_P,"with",100*a,"% confidence interval: [",CI[1],",",CI[2],"]\n")
}

#### Exe 6 ####
rm(list=ls())
set.seed(1)
n = 10000
Y=0
Y2 = 0
for(i in 1:n/2){
  U = runif(5)
  X = -log(U)
  X2 = -log(1-U)
  tmp2 = as.integer(((1:5)%*%X)>=21.6)
  tmp = (as.integer(((1:5)%*%X)>=21.6)+as.integer(((1:5)%*%X2)>=21.6))/2
  Y = Y + tmp
  Y2 = Y2 + tmp^2
}
Y/n
(Y2-Y^2/n)/(n-1)


#### Exe 8 ####
rm(list=ls())
set.seed(1)

a=0.05
z_a = qnorm(1-a/2)
done = 0
p = 50
n = 1
Y = 0
Y2 = 0
theta_hat = 0
var_hat = 0
while(!done){
  Z = rnorm(1)
  Y = Z^3*(exp(Z)-exp(-Z))/2
  Y2 = Y2 + Y^2
  theta_hat = theta_hat + (Y-theta_hat)/n
  if (n>1)
    var_hat = (Y2 - n*theta_hat^2)/(n-1)
  if((n >= p) & (2*sqrt(var_hat/n)*z_a <= 0.1))
    done = 1
  else
    n = n + 1
}
theta_hat
var_hat
c(theta_hat-z_a*sqrt(var_hat/n),theta_hat+z_a*sqrt(var_hat/n))


#### Exe 5 ####
library(Matrix)
library(matrixStats)
rm(list=ls())
set.seed(1)
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
