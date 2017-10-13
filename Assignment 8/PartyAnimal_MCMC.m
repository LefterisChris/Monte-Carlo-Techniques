function [] = PartyAnimal_MCMC()   % Code written by Suraj Keshri
% Party Animal: Exercise 3.1 from Barber

% x are the realizations (P D U H A), in this order
cur_x = [0 0 0 1 1];
cur_lik = PA_prob(cur_x);

Nsteps = 100000;

% Now run Metropolis-Hastings
mh_iter = zeros(Nsteps,1);
for n=1:Nsteps
    
    prp_x = cur_x;
    idx = randi(3);
    prp_x(idx) = ~prp_x(idx);
    prp_lik = PA_prob(prp_x);
    
    if rand < prp_lik/cur_lik
        cur_x = prp_x;
        cur_lik = prp_lik;
    end
    
    mh_iter(n) = cur_x(1);
end

P_avg_MH = sum(mh_iter)/Nsteps;

% Now run Gibbs
g_iters = zeros(Nsteps,1);
cur_x = [0 0 0 1 1];
for n=1:Nsteps
    for m=1:3
        p0 = PA_prob([cur_x(1:(m-1)) 0 cur_x((m+1):5)]);
        p1 = PA_prob([cur_x(1:(m-1)) 1 cur_x((m+1):5)]);
        if rand < p0/(p0+p1);
            cur_x(m) = 0;
        else
            cur_x(m) = 1;
        end
    end
    g_iters(n) = cur_x(1);
end

P_avg_G = sum(g_iters)/Nsteps;

% Now plot convergence of the two algorithm
figure;
subplot(2,1,1);
hold on;
plot(cumsum(mh_iter)./(1:Nsteps)');
ylim([0.5 0.7])
plot(0.609659561*ones(Nsteps,1),'k--')
grid on
hold off;
xlabel('Iteration');
legend('MH estimate','True probability value');

subplot(2,1,2);
hold on;
plot(cumsum(g_iters)./(1:Nsteps)');
plot(0.609659561*ones(Nsteps,1),'k--')  % The true answer
ylim([0.5 0.7])
grid on
hold off;
xlabel('Iteration');
legend('Gibbs estimate','True probability value');
end

%*************************************************************************
function p = PA_prob(x)

% x are the realizations (P D U H A), in this order

p = 0.2*(x(1)==1)+0.8*(x(1)==0);      % P
p = p*(0.4*(x(2)==1)+0.6*(x(2)==0));  % D
p = p*( (x(1)==1 && x(2)==1)*(0.999*(x(3)==1)+0.001*(x(3)==0)) + ...  
        (x(1)==1 && x(2)==0)*(0.9*(x(3)==1)+0.1*(x(3)==0)) + ... 
        (x(1)==0 && x(2)==1)*(0.9*(x(3)==1)+0.1*(x(3)==0)) + ...
        (x(1)==0 && x(2)==0)*(0.01*(x(3)==1)+0.99*(x(3)==0))); % U|P,D
p = p*( (x(1)==1)*(0.9*(x(4)==1)+0.1*(x(4)==0)) + ...
        (x(1)==0)*(0.2*(x(4)==1)+0.8*(x(4)==0))); % H|P
p = p*( (x(3)==1)*(0.95*(x(5)==1)+0.05*(x(5)==0)) + ...
        (x(3)==0)*(0.5*(x(5)==1)+0.5*(x(5)==0))); % A|U
end
