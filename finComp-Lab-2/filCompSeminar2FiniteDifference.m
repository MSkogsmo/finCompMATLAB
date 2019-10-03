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
M = 500; % Number of values we want to calculate, higher M gives higher precision
N = 10000; %Timesteps
delta_t = T/M;

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
%Här vet vi inte om det är så att r(t)-q(t) är bara r, var kommar gamma in?
beta_j_k = r*S_j*(delta_tao/(2*S_delta));

%l_j_k
l_j_k = alfa_j_k-beta_j_k;

%d_j_k
d_j_k = 1-r*delta_tao-2*alfa_j_k;

%u_j_k
u_j_k = alfa_j_k + beta_j_k;

%Setting tridiagonal matrix A
A1 = diag(d_j_k(2:M), 0);
A2 = diag(l_j_k(3:M), -1);
A3 = diag(u_j_k(2:M-1), 1);
A = A1+A2+A3;
%Values from outside the range exists in the matrix

%Calculating value function v_n and taking away negatives
v_n = (S_j - K)';
v_n(v_n < 0) = 0;


for i = N:-1:0
    v_n_temp = A*v_n(2:M); %Lägga till första värdet
    v_n_temp(M-1) = v_n_temp(M-1) + v_n(M+1)*u_j_k(M);
    v_n = [0; v_n_temp; u_j_k(M+1)*(S_max-K*exp(-r*(T-i*delta_t)))]; 
    
    %Använda s-k exp enligt slides, v_n(M+1) ska bytas ut mot randvärdet s-k exp från slisen
    %Kolla bok s.117-119, tror 119 mer han pekade på i figur 4.14 där vi
    %adderar ska (s-k)+exp(joxxx) in som sista värde istället för v_n(M+1) (den senare),
    %detta ska ske tiigare, inte som jag förstår det i läget där vi gör v_n
    %= (S_j - K)'; nu
    %Randvärdet
end
  
index_constant = S_0/S_max; %The relative place in the vector where we find the value
index = index_constant*M + 1; %The index for where S_0 exists
idx1 = ceil(index); %Upper value
idx2 = floor(index); %Lower vaue
interpolated_val = (v_n(idx2)+v_n(idx1))/2 %Extracting mean

plot(S_j,v_n)
%Bör plotta för alla värden för S_j
%Brde increasa

%Sänker N och M ger att explicita ej funkar, finns att behöver vara
%<>konstant
%Stabilitet slides och bok, vad händer om man minskar eller höjer värde
%utanför det som bör vara N och M

%Introduce grid for time and space

%Denote the value function to time and space

%Discriticize in forward or backward, 
%Implicit or explicit (stability value has not to do with accuracy, make sure not to big timesteps)
%delta_t <= c*power(delta_s,2) 
%v_delta = A*v_delta %Explicit time-stepping

%Insert in black-scholes algorithm

%Fix delta_t and delta_s to see what happens whan varying the other values