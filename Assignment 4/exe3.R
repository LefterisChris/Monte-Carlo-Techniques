#### exe 3 ####
set.seed(1)
N = 1000
lambda = 5
delta = 1/100
c = 600
A = 10000

S = rep(0,N)
t_A = rep(0,N)
for(i in 1:N){
  while(S[i] < A){
    U = runif(1)
    X_neg = log((1+delta*c/lambda)*U)/delta  # composition sampling
    if(X_neg < 0)
      X = X_neg
    else{
      X_pos = -log((1+lambda/delta/c)*(1-U))*c/lambda
      if(X_pos < 0 ){
        cat("Something went wrong")
        stop(1)
      }
      X = X_pos
    }
    S[i] = S[i] + X
    t_A[i] = t_A[i] +1
  }
}
Y = exp(-(delta-lambda/c)*S)
theta = mean(Y)
var_theta = var(Y)
a = 1-0.95
z = qnorm(1-a/2)
CI = c(theta - z*sqrt(var_theta/N),theta + z*sqrt(var_theta/N))
cat("Theta:",theta,"\nvariance:",var_theta,"\n95% confidence interval:",CI)
