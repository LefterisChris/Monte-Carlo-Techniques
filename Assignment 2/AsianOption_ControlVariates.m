% Estimation of c
n = 1000; m = 11; r=.05; sigma=.25; T=1; S0=100;    z_alpha = 2.575; 
K = [90,100,110,120]; P = zeros(n, 1);
mean=(r-(sigma^2)/2)*T/m; std =sigma*sqrt(T/m); S = zeros(m+1,1);

for i=1:n
    e_dX=exp(mean+std*randn(m,1));
    S(1)=S0;
    for j=1:m
        S(j+1)=S(j)*e_dX(j);
    end
    P(i) = sum(S(2:m+1))/m;
end
    
for j = 1 : max(size(K))
    Q = (P > K(j)).*(P - K(j));
    cov_est = cov([P Q]);
    c(j) = -cov_est(1,2)/cov_est(1,1);
end
disp('c(i)='),disp(c);

% Estimation of the option price
n = 10000; P = zeros(n, max(size(K)));

EX = 0;
for j = 1 : m
    EX = EX + exp(r*j*T/m);
end
EX = S0*EX/m;

for i=1:n
    e_dX=exp(mean+std*randn(m,1));
    S(1)=S0;
    for j=1:m
        S(j+1)=S(j)*e_dX(j);
    end
    avg = sum(S(2:m+1))/m;
    P(i,:) = (avg>K).*(avg-K) + c*(avg - EX);
    Naive_P(i,:) = (avg>K).*(avg-K);
end

for j = 1 : max(size(K))
    MU(j) = sum(P(:,j))/n;
    STD(j) = sqrt(sum((P(:,j)-MU(j)).^2)/(n-1));
    
    Naive_MU(j) = sum(Naive_P(:,j))/n;
    Naive_STD(j) = sqrt(sum((Naive_P(:,j)-Naive_MU(j)).^2)/(n-1));
end

disp('PRICE using control variate ='), disp(MU);
disp('PRICE using control variate ='), disp(Naive_MU);

CI = [MU - z_alpha*STD/sqrt(n);  MU + z_alpha*STD/sqrt(n)];
Naive_CI = [Naive_MU - z_alpha*Naive_STD/sqrt(n);  Naive_MU + z_alpha*Naive_STD/sqrt(n)];
disp('CI 99% using control variate ='), disp(CI);
disp('Naive CI 99% ='), disp(Naive_CI);