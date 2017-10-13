N = 100000;
tmp=0;
for i=1:N
    U = rand(1);
    tmp = tmp +U*exp(U)/N;
end

tmp-0.5*(exp(1)-1)
%%
U = rand(1,100000);
% cov(U,exp(U))
mean((exp(U)-exp(1)+1).*(U-0.5))
%%
clear all
clc

N = 10000;
M=100;
tmp = zeros(1,N);
for j=1:N
    cards = randperm(M);
    for i=1:M
        if i == cards(i)
            tmp(j) = tmp(j)+1;
        end
    end
end
avg=mean(tmp)
variance=var(tmp)


