%Fixed constants
K = 15;
S_0 = 14;
r = 0.1;
T = 0.5;
gamma = 1;
sigma = 0.25;

disc_errors_EM = [];
disc_errors_RK = [];

%Parameters
N = 10^(6);

%For plotting
sampled_prices = zeros(N,1);
disc_errors = [];
analytical = 0.8670;

for M=100:100:1500 %For loop for simulation variating M
    M %print M
    for i = 1:N %For loop for number of samples
        delta_t = T/M; %Rewriting delta_t
        S_i_1 = S_0; %Fixed variable for new calculations
        S_i_1_RK = S_i_1;
        S_i_1_EM = S_i_1;
        for j = 0:delta_t:T %For loop for simulating time-step
            %Runge-Kutta
            W_delta = randn*sqrt(delta_t);
            X_skatt = S_i_1_RK+r*S_i_1_RK*delta_t+sigma*S_i_1_RK*sqrt(delta_t);
            X_part1 = 1/(2*sqrt(delta_t));
            X_part2 = (sigma*X_skatt-sigma*S_i_1_RK)*(power(W_delta,2)-delta_t);
            S_i_1_RK = S_i_1_RK+r*S_i_1_RK*delta_t+sigma*S_i_1_RK*W_delta+(X_part1*X_part2);
            
            %Euler Maruyama
            S_i_1_EM = S_i_1_EM + r*delta_t * S_i_1_EM + sigma*power(S_i_1_EM,gamma)*(randn)*sqrt(delta_t);
        end
        sampled_prices_EM(i) = (max(S_i_1_EM - K, 0));  
        sampled_prices_RK(i) = (max(S_i_1_RK - K, 0));  
    end
    value_EM = mean(sampled_prices_EM)*exp(-r*T)
    value_RK = mean(sampled_prices_RK)*exp(-r*T)
    
    disc_errors_EM = [disc_errors_EM; abs(value_EM-analytical)];
    disc_errors_RK = [disc_errors_RK; abs(value_RK-analytical)];
end
loglog(disc_errors_EM)
hold on
loglog(disc_errors_RK)
hold off

