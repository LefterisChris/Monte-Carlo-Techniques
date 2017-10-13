clc;
n = 10000;
disp('True values')
theta = exp(-20); % Exact value
TrueVarEstimator = theta*(1-theta)/n; % Variance of estimator
disp('Theta = '), disp(theta)
disp('TrueVarEstimator = '), disp(TrueVarEstimator)

U = rand(n,1);
% Pure Monte Carlo simulation for rare event
disp('Naive Monte-Carlo')
X = (-log(U)) > 20;
NaiveEst = sum(X)/n; 
Naive_EstVar = sum((X-NaiveEst).^2)/(n-1);
disp('NaiveEst = '), disp(NaiveEst)
disp('Naive_EstVar = '), disp(Naive_EstVar)

%Importance sampling for rare event
disp('Importance Sampling for Rare Event')
t = 19/20;
Z = -log(U)/(1-t);
X = (Z > 20)./((1-t)*exp(t*Z));

ImpEst = sum(X)/n; 
Imp_EstVar = sum((X-ImpEst).^2)/(n-1);
disp('ImpEst = '), disp(ImpEst)
disp('Imp_EstVar = '), disp(Imp_EstVar)