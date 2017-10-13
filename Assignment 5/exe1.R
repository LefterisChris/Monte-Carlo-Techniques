#### Exe 1 ####
rm(list=ls())
# set.seed(1)

Blackscholes <- function(S,K,Dt,r,sigma) {
  d1 = (log(S/K)+(r+sigma^2/2)*Dt)/(sigma*sqrt(Dt))
  d2 = d1 - sigma*sqrt(Dt)
  return(S*pnorm(d1)-K*exp(-r*Dt)*pnorm(d2))
}


S0 = 100
K=100
Maturity = 0.5
r = 0.01
sigma = 0.4

BS_price = Blackscholes(S0,K,Maturity,r,sigma)

m = 10^seq(1,3,0.25)
p=length(m)
h = Maturity/m
mae1 = rep(0,p)
mae2 = rep(0,p)
# sample paths
N = 100000

for(s in 1:p){
  C0 = 0
  C0_2 = 0
  for(i in 1:N){
    X_hat = S0
    X_hat2 = S0
    Z = rnorm(m[s])
    for(j in 1:m[s]){
      # Euler scheme
      X_hat = X_hat + r*X_hat*h[s] + sigma*X_hat*sqrt(h[s])*Z[j]
      if(j %% 2 == 0){
        # Euler with Richardson extrapolation scheme
        X_hat2 = X_hat2 + r*X_hat2*2*h[s] + sigma*X_hat2*sqrt(h[s])*(Z[j-1]+Z[j])
      }
    }
    tmp1 = exp(-r*Maturity)*max(X_hat-K,0)
    tmp2 = exp(-r*Maturity)*max(X_hat2-K,0) 
    C0 = C0 + tmp1/N
    C0_2 = C0_2 + (2*tmp1-tmp2)/N
  }
  cat("C0 = ",C0,"C0 (2) = ",C0_2,"\n")
  mae1[s] = abs(BS_price-C0)
  mae2[s] = abs(BS_price-C0_2)
}

plot(m,mae1,log="xy",type="l",lty=2,ylab = "Absolute error",main = "Euler vs Euler /w Richardson scheme")
lines(m,mae2,log = "xy",type = "l",lty=3)
grid(equilogs = F)
legend("bottomleft",lty = c(2,3),legend = c("Euler","Euler /w Richardson"))
