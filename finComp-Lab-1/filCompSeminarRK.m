%Fixed constants
K = 15;
S_0 = 14;
r = 0.1;
T = 0.5;
gamma = 1;
sigma = 0.25;

%Parameters
M = 200;
N = 10^(7);
delta_t = T/M;

%For plotting
sampled_prices = zeros(N,1);
disc_errors = [];
analytical = 0.8670;

for gamma=0:20:200 %For loop for simulation
    for i = 1:N %For loop for number of samples
        delta_t = T/n; %Rewriting delta_t
        S_i_1 = S_0; %Fixed variable for new calculations
        for j = 0:delta_t:T %For loop for simulating time-step
            %Runge-Kutta
            W_delta = randn*sqrt(delta_t);
            X_skatt = S_i_1+r*S_i_1*delta_t+sigma*S_i_1*sqrt(delta_t);
            X_part1 = 1/(2*sqrt(delta_t));
            X_part2 = (sigma*X_skatt-sigma*S_i_1)*(power(W_delta,2)-delta_t);
            S_i_1 = S_i_1+r*S_i_1*delta_t+sigma*S_i_1*W_delta+(X_part1*X_part2);
        end
        sampled_prices(i) = (max(S_i_1 - K, 0));    
    end
    value = mean(sampled_prices)*exp(-r*T)
    disc_errors = [disc_errors; abs(value-analytical)];
end
loglog(disc_errors)

