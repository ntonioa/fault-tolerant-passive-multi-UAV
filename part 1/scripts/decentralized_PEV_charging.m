clear, clc, close all;

%% Parameters
% Simulation
N = 24; % Prediction horizon
m = 500; % Number of PEVs
T = 1/3; % Time step duration
k_max = 100; % Maximum number of iterations

% Global
P_max_1 = 500*ones(1, 10);
P_max_2 = 200*ones(1, 4); 
P_max = [P_max_1, P_max_2, P_max_1]; % Concatenated maximum power for each time step
P_ref = [290*ones(1, N-1), 0]; % Power reference for each time step (except the last one)
tol = 25; % Tolerance for power reference deviation

% Local
F = N*ones(m, 1); % Final step for each PEV
x_init = zeros(m, 1); % Initial state of charge for each PEV
pevs(1:m, 1) = PevMpc; % Array of PEV objects
for p = 1:m
    x_max = 8*(1+rand); % Maximum state of charge for each PEV
    x_init(p) = (0.2+0.3*rand)*x_max; % Initial state of charge for each PEV
    x_ref = (0.55+0.25*rand)*x_max; % State of charge reference for each PEV
    eta_ch = 0.925+0.06*rand; % Charging efficiency for each PEV
    eta_dis = 0.925+0.06*rand; % Discharging efficiency for each PEV
    xi = 0.3*rand(1, N); % Random perturbation for each PEV

    pevs(p, 1) = PevMpc(N, T, x_max, 1, x_ref, 5, 1.3, 0, 0, eta_ch, 1, xi); % Create a PevMpc object for each PEV
end

%% Variables
P = zeros(m, N, k_max+1); % Power profile for each PEV at each time step and iteration

rho = zeros(m, N, k_max); % Variable rho for each PEV at each time step and iteration
rho_agg = zeros(N, k_max); % Aggregated variable rho for all PEVs at each time step and iteration
lambda = zeros(N, k_max); % Multiplier lambda for power constraint at each time step and iteration
mu = zeros(N, k_max); % Multiplier mu for power reference at each time step and iteration
nu = zeros(N, k_max); % Multiplier nu for power reference at each time step and iteration
alpha = 0.001./((2:k_max)); % Step size for updating multipliers

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
    mu_curr = mu(:, k); % Current value of mu
    ni_curr = nu(:, k); % Current value of nu
    time = 0; % Total time for parallel computing
    parfor p = 1:m
        tic; % Start timer
        pevs(p) = pevMpcIter(pevs(p), x_init(p), F(p), lambda_curr, mu_curr, ni_curr); % Perform MPC iteration for each PEV
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

    rho(:, :, k+1) = reshape([pevs.s_up]-[pevs.s_down], m, []); % Update variable rho
    rho_agg(:, k+1) = N*max(rho(:, :, k+1), [], 1); % Update aggregated variable rho
    lambda(:, k+1) = max(zeros(N, 1), lambda(:, k)+alpha(k)*(sum(P(:, :, k+1), 1)'-P_max'+rho_agg(:, k+1))); % Update multiplier lambda
    mu(:, k+1) = max(zeros(N, 1), mu(:, k)+alpha(k)*(sum(P(:, :, k+1), 1)'-P_ref')); % Update multiplier mu
    nu(:, k+1) = max(zeros(N, 1), nu(:, k)-alpha(k)*(sum(P(:, :, k+1), 1)'-P_ref')); % Update multiplier nu

    % Plot aggregated power step-by-step
    stairs(0:T:(N-1)*T, P_max, 'r'); % Plot maximum power limit
    hold on;
    stairs(0:T:(N-1)*T, P_ref, '--'); % Plot power reference
    stairs(0:T:(N-1)*T, sum(P(:, :, k+1), 1)); % Plot current aggregated power
    drawnow;
    hold off;
end

close;
%% Plots
% Aggregated power
P_agg = sum(P(:, :, k+1), 1); % Aggregated power profile
figure, hold on, grid on;
stairs(0:T:(N-1)*T, P_agg, 'LineWidth', 1.5); % Plot aggregated power
stairs(0:T:(N-1)*T, P_max, '--', 'LineWidth', 1.5); % Plot maximum power limit
stairs(0:T:(N-1)*T, P_ref, ':', 'LineWidth', 2); % Plot power reference
title('\textbf{Aggregated power profile}', 'Interpreter', 'LaTeX'); 
xlabel('Time [$h$]', 'Interpreter', 'LaTeX'); 
ylabel('Power [$kW$]', 'Interpreter', 'LaTeX'); 
xlim([0 (N-1)*T-0.01]), ylim([0 520]);
legend('Aggregated power', 'Maximum power', 'Reference power', 'Location', 'south', 'Interpreter', 'LaTeX');

% Global constraint violation
P_max_viol = squeeze(max(sum(P(:, :, 2:k+1), 1)-P_max)); % Maximum violation of the maximum aggregated power constraint
P_ref_gap = squeeze(max(abs(sum(P(:, :, 2:k+1), 1)-min(P_ref, P_max)))); % Maximum deviation from the adjusted aggregated power reference
figure;

subplot(2, 1, 1);
plot(0:k-1, P_max_viol, 'LineWidth', 1.5); % Plot maximum power constraint violation
yline(0, '--r', 'LineWidth', 1.5); % Plot threshold
title('\textbf{Maximum aggregated power constraint violation}', 'Interpreter', 'LaTeX'); 
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('Violation [$kW$]', 'Interpreter', 'LaTeX'); 
xlim([0 k-1]), ylim([min(P_max_viol)-20 max(P_max_viol)+20]);
legend('Maximum constraint violation', 'Threshold', 'Location', 'northeast', 'Interpreter', 'LaTeX');
grid on;

subplot(2, 1, 2);
plot(0:k-1, P_ref_gap, 'LineWidth', 1.5); % Plot power reference gap
yline(tol, '--r', 'LineWidth', 1.5); % Plot threshold
title('\textbf{Maximum aggregated power reference gap}', 'Interpreter', 'LaTeX'); 
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('Gap [$kW$]', 'Interpreter', 'LaTeX'); 
xlim([0 k-1]), ylim([min(P_ref_gap)-20 max(P_ref_gap)+20]);
legend('Maximum reference gap', 'Threshold', 'Location', 'northeast', 'Interpreter', 'LaTeX');
grid on;

% State
figure, hold on, grid on;
rand_p = randi(m, 1, 10); % Randomly select 10 PEVs
for p = rand_p
    plot(0:T:N*T, pevs(p).sol.x, 'LineWidth', 1.5); % Plot state of charge for each selected PEV
end
title('\textbf{PEV state of charge}', 'Interpreter', 'LaTeX'); 
xlabel('Time [$h$]', 'Interpreter', 'LaTeX'); 
ylabel('Charge [$kWh$]', 'Interpreter', 'LaTeX'); 
xlim([0 N*T]), ylim([0 13]);

% Norm of rho
figure, hold on, grid on;
for p = rand_p
    plot(0:k-1, vecnorm(squeeze(rho(p, :, :))), 'LineWidth', 1.5); % Plot norm of variable rho for each selected PEV
end
title('\textbf{Norm of the variable $\rho$}', 'Interpreter', 'LaTeX'); 
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('$||\rho||$', 'Interpreter', 'LaTeX'); 
xlim([0 k-1]);

% Lambda
figure;
plot(0:k-1, lambda(1:N-1, :)', 'LineWidth', 1.5); % Plot multipliers lambda
title('\textbf{Multiplier $\lambda$}', 'Interpreter', 'LaTeX'); 
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('$\lambda$', 'Interpreter', 'LaTeX'); 
xlim([0 k-1]);
grid on;

% Mu
figure;
plot(0:k-1, mu(1:N-1, :)', 'LineWidth', 1.5); % Plot multipliers mu
title('\textbf{Multiplier $\mu$}', 'Interpreter', 'LaTeX'); 
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('$\mu$', 'Interpreter', 'LaTeX'); 
xlim([0 k-1]);
grid on;

% Nu
figure;
plot(0:k-1, nu(1:N-1, :)', 'LineWidth', 1.5); % Plot multipliers nu
title('\textbf{Multiplier $\nu$}', 'Interpreter', 'LaTeX'); 
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('$\nu$', 'Interpreter', 'LaTeX'); 
xlim([0 k-1]);
grid on;
