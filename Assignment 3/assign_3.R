N = 100000

U = runif(N)
Y = - 2*log(U)
X = rep(0,N)
Y_v = (Y >= 0 & Y<=4)
X[Y_v] = exp(-(4-Y[Y_v]))
X[!Y_v] = rep(1,sum(Y_v==FALSE))

mean(X)
var(X)

U2 = runif(N)
X2 = -log(U2)
Y2 = rep(0,N)
X_v = ( X2>=0 &  X2<= 4)
Y2[X_v] = exp(-0.5*(4-X2[X_v]))
Y2[!X_v] = rep(1,sum(X_v==FALSE))

mean(Y2)
var(Y2)

#### Exe 2 ####
rm(list=ls())
BlackScholes <- function(S0,K,s,r,Dt) {
  d1 = (log(S0/K)+(r+s^2/2)*Dt)/(s*sqrt(Dt))
  d2 = d1 - s*sqrt(Dt)
  return(S0*pnorm(d1)-K*exp(-r*Dt)*pnorm(d2))
}

library(matrixStats)
set.seed(1)
S0 = 100
Maturity = 1
r = 0.05
s=0.25
L=105
K1 = 110
K2 = 120

N = 10000

# usual simulation
St = matrix(0,nrow=N,ncol=3)
St[,1] = S0
St[,2] = St[,1]*exp((r-s^2/2)*Maturity/2+s*sqrt(Maturity/2)*rnorm(N))
St[,3] = St[,2]*exp((r-s^2/2)*Maturity/2+s*sqrt(Maturity/2)*rnorm(N))

K = rep(K1,N)
Indicator = (St[,2] > L) 
K[Indicator] = K2

payoff = rowMaxs(matrix(c(St[,3]-K,rep(0,N)),nrow=N))
C0 = mean(exp(-r*Maturity)*payoff)
# 95% confidence interval
s_var = sd(exp(-r*Maturity)*payoff)
z = qnorm(1-0.05/2)
CI = c(C0-z*s_var/sqrt(N),C0+z*s_var/sqrt(N))
C0
s_var
CI

# conditional expectation
BSprices = exp(-r*Maturity/2)*BlackScholes(St[,2],K,s,r,Maturity/2)
BS_pr = mean(BSprices)
s_var_BS = sd(BSprices)
CI_BS = c(BS_pr-z*s_var_BS/sqrt(N),BS_pr+z*s_var_BS/sqrt(N))
BS_pr
s_var_BS
CI_BS

# using antithetic vars as well
set.seed(1)
St1 = matrix(0,nrow=N,ncol=2)
St1[,1] = S0
Z = rnorm(N)
St1[,2] = St1[,1]*exp((r-s^2/2)*Maturity/2)*(exp(s*sqrt(Maturity/2)*Z)+exp(s*sqrt(Maturity/2)*(-Z)))/2

K_av = rep(K1,N)
Indicator2 = (St1[,2] > L) 
K_av[Indicator2] = K2
BSprices2 = exp(-r*Maturity/2)*BlackScholes(St1[,2],K_av,s,r,Maturity/2)
BS_pr2 = mean(BSprices2)
s_var_BS2 = sd(BSprices2)
CI_BS2 = c(BS_pr2-z*s_var_BS2/sqrt(N),BS_pr2+z*s_var_BS2/sqrt(N))
BS_pr2
s_var_BS2
CI_BS2

#### exe4 ####
rm(list=ls())
N=10000

U = runif(N)
X = -log(U)
I = as.integer(X >20)
mean(I)
var(I)

U1 = runif(N)
X2 = -21*log(U)
I2 = as.integer(X2 > 20)
Y = I2 * exp(-X2)*21/exp(-X2/21)
mean(Y)
var(Y)
