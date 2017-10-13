clc;
% Initial input data
n = 10000;     m = 2;  r = .05;    sigma = .25;    T = 1;  S0 = 100;
dT = T/m; % Time increment between trading points
L = 105;    K1 = 110;   K2 = 120;    z_alpha = 1.96;
drift = (r-(sigma^2)/2)*dT; % Mean of GBM during trading
Vol = sigma*sqrt(dT); % Standard deviation of GBM during trading

S1 = S0*exp(drift + Vol*randn(1,n));
K = (S1 <= L)*K1+(S1 > L)*K2;

% Price computations for standard Monte-Carlo
disp('Standard Monte-Carlo')
S2 = S1.*exp(drift + Vol*randn(1,n));
DiscPayoff = exp(-r*T)*max(0, S2 - K);
PriceEst = mean(DiscPayoff);
disp('PriceEst ='), disp(PriceEst);
Std_DiscPayoff= std(DiscPayoff)/sqrt(n);
Standard_CI = [PriceEst - z_alpha*Std_DiscPayoff   PriceEst + z_alpha*Std_DiscPayoff];
disp('Standard CI ='), disp(Standard_CI);

% Price computations for conditional Monte-Carlo
disp('Conditional Monte-Carlo')
d2 = (log(S1./K) + drift)/(sigma*sqrt(dT));
d1 = d2 + sigma*sqrt(dT);
C1 = S1.*normcdf(d1) - exp(-r*dT)*K.*normcdf(d2);
C0 = exp(-r*dT)*C1;
Cond_PriceEst = sum(C0)/n;
disp('Cond_PriceEst ='), disp(Cond_PriceEst);
Cond_Std = std(C0)/sqrt(n);
Conditional_CI = [Cond_PriceEst - z_alpha*Cond_Std Cond_PriceEst + z_alpha*Cond_Std];
disp('Conditional CI ='), disp(Conditional_CI);

