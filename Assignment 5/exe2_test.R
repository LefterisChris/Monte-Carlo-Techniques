Maturity=1
S0 = K = 100
r=.05
V0=.04
k=1.2
theta = 0.04
cor = -0.5
sigma = 0.3

# step length
h = 0.01
# time steps
m = floor(Maturity/h)
# number of paths
N = 100000

C0 = 0.0
C0_2 = 0.0
for(i in 1:N){
  S_hat = S0
  V_hat = V0
  S_hat2 = S0
  V_hat2 = V0
  Z1 = rnorm(m)
  Z = rnorm(m)
  Z2 = cor*Z1+sqrt(1-cor^2)*Z
  for(j in 1:m){
    S_hat = S_hat + r*S_hat*h+sqrt(V_hat*h)*S_hat*Z1[j]
    V_hat = V_hat + k*(theta-V_hat)*h + sigma*sqrt(V_hat*h)*Z2[j]
    # V_hat = abs(V_hat)
    if(V_hat < 0){
      V_hat = 0
    }
    # if(j %% 2==0){
    #   S_hat2 = S_hat2 + r*S_hat2*2*h+sqrt(V_hat2*h)*S_hat2*(Z1[j/2]+Z1[j/2+1])
    #   V_hat2 = V_hat2 + k*(theta-V_hat2)*2*h + sigma*sqrt(V_hat2*h)*(Z2[j/2]+Z2[j/2+1])
    #   V_hat2 = abs(V_hat2)
    # }
  }
  for(j in 1:(m/2)){
      S_hat2 = S_hat2 + r*S_hat2*2*h+sqrt(V_hat2*h)*S_hat2*(Z1[2*j-1]+Z1[2*j])
      V_hat2 = V_hat2 + k*(theta-V_hat2)*2*h + sigma*sqrt(V_hat2*h)*(Z2[2*j-1]+Z2[2*j])
      # V_hat2 = abs(V_hat2)
      if(V_hat2 < 0){
        V_hat2 = 0
      }
  }
  tmp = exp(-r*T)*max(S_hat-K,0)
  tmp2 = exp(-r*T)*max(S_hat2-K,0)
  C0 = C0 + tmp/N
  C0_2 = C0_2 + (2*tmp-tmp2)/N
}
C0
C0_2
