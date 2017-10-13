#### exe 5 ####
rm(list=ls())
# set.seed(1)
r = 0.05
sigma = 0.25
Maturity = 1
S0 = 100
K = 90
m = 6

a = exp((r-sigma^2/2)*Maturity/m)
b = sigma*sqrt(Maturity/m)

# construct mu - factor for Z
mu = S0*b/m*rev(cumsum(rev(a^(1:m))))
mu = mu/norm(x = mu,type = "2")

N = 1000
M = 10
p=1/M

xs = qnorm(0:M/M)

price_st = 0
var_st = 0
for(i in 1:M){
  tmp1 = 0
  tmp2 = 0
  for(j in 1:N){
    # generate W
    U = runif(1)
    phi_w1 = pnorm(xs[i])
    phi_w2 = pnorm(xs[i+1])
    W = qnorm(U*(phi_w2-phi_w1)+phi_w1)
    Z = W*mu + (diag(m)-mu%*%t(mu))%*%rnorm(m)
    Z = as.vector(Z)
    payoff = exp(-r*Maturity)*max(S0/m*sum(a^(1:m)*exp(b*cumsum(Z)))-K,0)
    tmp1 = tmp1 + payoff
    tmp2 = tmp2 + payoff^2
  }
  theta_i = tmp1/N
  var_i = (tmp2 - N*theta_i^2)/(N-1)
  price_st = price_st + p*theta_i
  var_st = var_st + p^2*var_i/N
}
a = 0.99
z = qnorm((1+a)/2)
CI = c(price_st-z*sqrt(var_st),price_st+z*sqrt(var_st))
cat("Uniform Stratification\n")
cat("Option price:",price_st,"\nvariance:",var_st,"\n99% confidence interval:",CI,"\n")

accu = diff(CI)

# Regular MC to price the option
cnt = 1
theta = 0
sample_var = 1000
sum_reg = 0
sum_sqr = 0
not_done = T
while(not_done){
  GBM = rep(0,m+1)
  GBM[1] = S0
  W = rnorm(m)
  for(i in 1:m){
    GBM[i+1] = GBM[i]*exp((r-sigma^2/2)*(Maturity/m)+sigma*sqrt(Maturity/m)*W[i])
  }
  payoff = exp(-r*Maturity)*max(sum(GBM[-1])/m-K,0)
  sum_sqr = sum_sqr + payoff^2
  sum_reg = sum_reg + payoff
  theta = sum_reg/cnt
  if(cnt > 1){
    sample_var = (sum_sqr - cnt*theta^2)/(cnt-1)
  }
  cnt = cnt + 1
  lb = theta - z*sqrt(sample_var/cnt)
  ub = theta + z*sqrt(sample_var/cnt)
  if((ub-lb <= accu) & (cnt > 10000))
    not_done = F
}
theta
sample_var
cat("Number of sample paths required in the reg MC",cnt,"\n")

