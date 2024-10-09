%#ok<*PFBNS>

clear, clc, close all;
tol = 5e-2;
rng(202405);

%% Parameters
N = 24;
T = 1/3;
m = 500;
k_max = 25;
l_max = 300;

P_max = [0.8*m*ones(5, 1); 0.4*m*ones(10, 1); 1.2*m*ones(9, 1)];
% P_max = [0.5*m; inf*ones(23, 1)];

F = N*ones(1, m);
x_init = zeros(1, m);
pevs(1, 1:m) = PevMpc;
for p = 1:m
    x_max = 8*(1+rand);
    x_init(p) = (0.2+0.3*rand)*x_max;
    x_ref = (0.55+0.25*rand)*x_max;
    eta_ch = 0.925+0.06*rand;
    eta_dis = 0.925+0.06*rand;
    xi = 0.3*rand(N, 1);
    
    pevs(p) = PevMpc(N, T, x_max, 1, x_ref, 5, 1.3, 0, 0, eta_ch, eta_dis, xi);
end

alpha = 0.01./((1:l_max));

%% Variables
P = zeros(N, m, k_max+1); % Power profile for each PEV at each time step and iteration

rho = zeros(N, m, k_max); % Variable rho for each PEV at each time step and iteration
rho_agg = zeros(N, k_max); % Aggregated variable rho for all PEVs at each time step and iteration
lambda = zeros(N, k_max); % Multiplier lambda for power constraint at each time step and iteration
%% Solution
figure;

k = 0; % Iteration counter
while(true)
    k = k+1; % Increment iteration counter
    disp(" ");
    disp("Iteration "+(k-1));
    
    % Parallel computing
    P_next = zeros(m, N); % Placeholder for next power profile for each PEV
    lambda_curr = lambda(:, k); % Current value of lambda
    time = 0; % Total time for parallel computing
    parfor p = 1:m
        tic; % Start timer
        pevs(p) = pevMpcIter(pevs(p), x_init(p), F(p), lambda_curr, mu_curr, nu_curr); % Perform MPC iteration for each PEV
        P_next(p, :) = pevs(p).sol.P; % Store the next power profile for each PEV
        time = time+toc; % Accumulate time for parallel computing
    end
    P(:, :, k+1) = P_next; % Update power profiles
    disp("Average iteration time (per PEV): "+(time/m)+" s");
    % ------------------

    P_max_viol = max(sum(P(:, :, k+1), 1)-P_max); % Maximum violation of the maximum aggregated power constraint
    P_ref_gap = max(abs(sum(P(:, :, k+1), 1)-min(P_ref, P_max))); % Maximum deviation from the adjusted aggregated power reference
    disp("Maximum violation of the maximum aggregated power constraint: "+P_max_viol+" kW");
    disp("Maximum deviation from the adjusted aggregated power reference: "+P_ref_gap+" kW");
    if  (P_max_viol <= 0 && P_ref_gap <= tol) || k == k_max 
        P = P(:, :, 1:k+1); % Trim the power profiles to the final iteration
        rho = rho(:, :, 1:k); % Trim the variable rho to the final iteration
        rho_agg = rho_agg(:, 1:k); % Trim the aggregated variable rho to the final iteration
        lambda = lambda(:, 1:k); % Trim the multipliers lambda to the final iteration
        mu = mu(:, 1:k); % Trim the multipliers mu to the final iteration
        nu = nu(:, 1:k); % Trim the multipliers nu to the final iteration
        break; % Exit the loop
    end

    lambda(:, k+1) = max(zeros(N, 1), lambda(:, k)+alpha(k)*(sum(P(:, :, k+1), 1)'-P_max'+rho_agg(:, k+1))); % Update multiplier lambda

    % Plot aggregated power step-by-step
    stairs(0:T:(N-1)*T, P_max, 'r'); % Plot maximum power limit
    hold on;
    stairs(0:T:(N-1)*T, P_ref, '--'); % Plot power reference
    stairs(0:T:(N-1)*T, sum(P(:, :, k+1), 1)); % Plot current aggregated power
    drawnow;
    hold off;
end

close;