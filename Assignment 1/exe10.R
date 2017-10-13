rm(list=ls())
# set.seed(1)
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
group <- round(seq(from=floor(n/m),to=n, length.out = m))

for(i in 1:length(K)){
  payoff <-0
  for(j in 1:P){
    priceT <- sum(GBM[j,group])/m
    payoff <- payoff + max(0,priceT-K[i])
  }
  cat("For strike K=",K[i],":",exp(-r*T)*payoff/P,"\n")
}


# Brownian motion
# N=10000
# M <- 10
# BM = matrix(0,nrow=N,ncol=M)
# BM[,1] <- rep(0,N)
# for(i in 1:(M-1)){
#   BM[,i+1] <- BM[,i]+sqrt(1/M)*rnorm(N)
# }

# CV <- cov(BM)
# C <- t(chol(CV[2:10,2:10]))
# 
# BM2 <- C%*%matrix(rnorm(M-1),ncol=1)
# 
# s2 <- matrix(0,nrow=M,ncol=M)
# for(i in 1:M)
#   s2[i,1:i]<-sqrt(1/M)
# 
# BM2 <- s2*rnorm(N)
# 
# 
# 
# 
# 
