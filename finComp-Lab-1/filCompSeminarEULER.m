%Fixed constants
K = 15;
S_0 = 14;
r = 0.1;
T = 0.5;
gamma = 1;
sigma = 0.25;

%Parameters
M = 200;
N = 10^(6);
delta_t = T/M;

%For plotting
sampled_prices = zeros(N,1);
disc_errors = [];
analytical = 0.8670;

for gamma=0.5:0.05:1 %For loop for simulation
    for i = 1:N %For loop for number of samples
        delta_t = T/N; %Rewriting delta_t
        S_i_1 = S_0; %Fixed variable for new calculations
        for j = 0:delta_t:T %For loop for simulating time-step
            %Euler?Maruyama
            S_i_1 = S_i_1 + r*delta_t * S_i_1 + sigma*power(S_i_1,gamma)*(randn)*sqrt(delta_t);
        end
        sampled_prices(i) = (max(S_i_1 - K, 0));    
    end
    value = mean(sampled_prices)*exp(-r*T)
    disc_errors = [disc_errors; value];
    
end

plot(disc_errors)

