%Implicit method

interpolated_val_array = [];
for M = 0:100:1000
    %Constants - Same as the other other seminar but with finance difference
    %Fixed constants
    K = 15;
    S_0 = 14;
    r = 0.1;
    T = 0.5;
    gamma = 1;
    sigma = 0.25;
    analytical = 0.8670; %Analytical value from bsect calc

    %Parameters, här är vi osäkra om hur det ligger till med hur matrixen
    %blir...
    N = 10000; %Timesteps
    delta_t = T/N;

    %Calculating formulas, setting up matricies

    %S_j
    S_min = 0; %Setting lower boundry
    S_max = 4*K; %Setting upper boundry
    S_delta = (S_max-S_min)/M;
    S_j = [0:M]*S_delta + S_min;

    %Tao
    delta_tao = T/N;
    Tao_k = [0:N]*delta_tao; %This is not used in the calculations

    %alfa
    alfa_j_k = (delta_tao/(power(S_delta,2)))*((power(sigma,2)*power(S_j,2))/2);

    %beta
    beta_j_k = r*S_j*(delta_tao/(2*S_delta));

    %l_j_k
    l_j_k = -alfa_j_k+beta_j_k;

    %d_j_k
    d_j_k = 1+r*delta_tao+2*alfa_j_k;

    %u_j_k
    u_j_k = -alfa_j_k - beta_j_k;

    %Setting tridiagonal matrix A
    A1 = diag(d_j_k(1:M+1), 0);
    A2 = diag(l_j_k(2:M+1), -1);
    A3 = diag(u_j_k(1:M), 1);
    A = A1+A2+A3; %Values from outside the range exists in the matrix
    A = inv(A);
    
    %Calculating value function v_n and taking away negatives
    v_n = (S_j - K)';
    v_n(v_n < 0) = 0;
    
    for i = N:-1:0
        v_n_temp = A*v_n;
        v_n_temp(M+1) = S_max-K*exp(-r*(T-i*delta_t));
        v_n = v_n_temp;
    end
    index_constant = S_0/S_max; %The relative place in the vector where we find the value
    index = index_constant*M + 1; %The index for where S_0 exists
    idx1 = ceil(index); %Upper value
    idx2 = floor(index); %Lower vaue
    interpolated_val = (v_n(idx2)+v_n(idx1))/2; %Extracting mean
    interpolated_val_array = [interpolated_val_array, interpolated_val];
    plot(S_j,v_n);
    hold on;
end

interpolated_val_array

plot(S_j,v_n)
