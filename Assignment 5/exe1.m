clc
clear all
close all
%% Exe 1 
S0 = 100;
K = 100;
Maturity = 0.5;
r = 0.01;
sigma = 0.4;

BS_price = BlackScholes(S0,K,Maturity,sigma,r);

m = round(10.^(1:0.25:3));
p = length(m);
h = Maturity./m;

mae1 = zeros(1,p);
mae2 = zeros(1,p);
% sample paths
N = 100000;

for s=1:p
  C0 = 0;
  C0_2 = 0;
  for i=1:N
    X_hat = S0;
    X_hat2 = S0;
    Z = randn(m(s));
    for j=1:m(s)
      % Euler scheme
      X_hat = X_hat + r*X_hat*h(s) + sigma*X_hat*sqrt(h(s))*Z(j);
    end
    for j=1:floor(m(s)/2)
      % Euler with Richardson extrapolation scheme
      X_hat2 = X_hat2 + r*X_hat2*2*h(s) + sigma*X_hat2*sqrt(h(s))*(Z(2*j-1)+Z(2*j));
    end
    tmp1 = exp(-r*Maturity)*max(X_hat-K,0);
    tmp2 = exp(-r*Maturity)*max((2*X_hat2-X_hat)-K,0);
    C0 = C0 + tmp1/N;
    C0_2 = C0_2 + tmp2/N;
  end
  disp(['C0= ',num2str(C0),' C0_2 = ',num2str(C0_2),'\n'])
  mae1(s) = abs(BS_price-C0);
  mae2(s) = abs(BS_price-C0_2);
end

