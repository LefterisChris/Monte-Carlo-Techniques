rm(list=ls())
lambda_t<-function(t){
  return((1+2*t+t^2)/100)
}

lambda = 121/100;
T=10

expenses <- 0
N <- 10000
for(i in 1:N){
  I = 0
  S <- c()
  U1 = runif(1)
  t = 0
  t = t - log(U1)/lambda
  while(t<T){
    U2 = runif(1)
    if(U2 <= lambda_t(t)/lambda){
      I <- I + 1
      S <- c(S,t)
      U3 <- runif(1)
      if(U3 <= 1/3){
        expenses <- expenses + 100
      }else if(U3 <= 2/3){
        expenses <- expenses + 500
      }else
        expenses <- expenses + 1000
    }
    U1 = runif(1)
    t <- t - log(U1)/lambda
  }
}
cat(expenses/N)