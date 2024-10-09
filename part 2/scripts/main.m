clear, clc, close all;

%% Parameters
% Simulation
N = 24; % Prediction horizon
T = 1/3; % Time step duration
m = 20; % Number of PEVs
k_max = 4; % Maximum number of iterations

% Global
P_max_1 = 20.0*ones(10, 1);
P_max_2 = 500.0*ones(4, 1);
P_max = [P_max_1; P_max_2; P_max_1]; % Concatenated maximum power for each time step
P_ref = [2.90*ones(N-1, 1); 0]; % Power reference for each time step (except the last one)
% tol = 25; % Tolerance for power reference deviation

% Local
F = N*ones(1, m); % Final step for each PEV
x_init = zeros(1, m); % Initial state of charge for each PEV
pevs(1, 1:m) = PevMpc; % Array of PEV objects
for p = 1:m
    x_max = 8*(1+rand); % Maximum state of charge for each PEV
    x_init(p) = (0.2+0.3*rand)*x_max; % Initial state of charge for each PEV
    x_ref = (0.55+0.25*rand)*x_max; % State of charge reference for each PEV
    eta_ch = 0.925+0.06*rand; % Charging efficiency for each PEV
    eta_dis = 0.925+0.06*rand; % Discharging efficiency for each PEV
    xi = 0.3*rand(N, 1); % Random perturbation for each PEV

    pevs(p) = PevMpc(N, T, x_max, 1, x_ref, 5, 1.3, 0, 0, eta_ch, 1, xi); % Create a PevMpc object for each PEV
end

%% Variables
[P_tilde_star, rho] = algorithm1(P_max, P_ref, pevs, N, m, k_max, F, x_init);

%%
P_agg = sum(P_tilde_star, 2); % Aggregated power profile
figure, hold on, grid on;
stairs(0:T:(N-1)*T, P_agg, 'LineWidth', 1.5); % Plot aggregated power
stairs(0:T:(N-1)*T, P_max, '--', 'LineWidth', 1.5); % Plot maximum power limit
stairs(0:T:(N-1)*T, P_ref, ':', 'LineWidth', 2); % Plot power reference
title('\textbf{Aggregated power profile}', 'Interpreter', 'LaTeX');
xlabel('Time [$h$]', 'Interpreter', 'LaTeX');
ylabel('Power [$kW$]', 'Interpreter', 'LaTeX');
xlim([0 (N-1)*T-0.01]), ylim([0 520]);
legend('Aggregated power', 'Maximum power', 'Reference power', 'Location', 'south', 'Interpreter', 'LaTeX');
