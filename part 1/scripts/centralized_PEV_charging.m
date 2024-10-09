clear, clc, close all;

%% Parameters
% Simulation
N = 24; % Prediction horizon
m = 500; % Number of PEVs
T = 1/3; % Time step duration

% Global
P_max_1 = 500*ones(1, 10);
P_max_2 = 200*ones(1, 4);
P_max = [P_max_1, P_max_2, P_max_1]; % Concatenated maximum power for each time step
P_ref = [290*ones(1, N-1), 0]; % Reference power for all time steps

% Local
F = N*ones(m, 1); % Final step for each PEV
P_ch_max = 5*ones(m, 1); % Maximum charging power for each PEV
P_ch_min = 1.3*ones(m, 1); % Minimum charging power for each PEV
P_dis_max = 0*ones(m, 1); % Maximum discharging power for each PEV
P_dis_min = 0*ones(m, 1); % Minimum discharging power for each PEV
x_max = 8*(ones(m, 1)+rand(m, 1)); % Maximum state of charge for each PEV
x_min = 1*ones(m, 1); % Minimum state of charge for each PEV
x_init = (0.2*ones(m, 1)+0.3*rand(m, 1)).*x_max; % Initial state of charge for each PEV
x_ref = (0.55*ones(m, 1)+0.25*rand(m, 1)).*x_max; % Reference state of charge for each PEV
eta_ch = 0.925*ones(m, 1)+0.06*rand(m, 1); % Charging efficiency for each PEV
eta_dis = 1*ones(m, 1); % Discharging efficiency for each PEV
xi = 0.3*rand(m, N); % Random perturbation for each PEV and time step

%% Solution
tic
sol = CenMpc(N, T, m, x_min, x_max, x_init, x_ref, F, P_ch_min, P_ch_max, P_dis_min, P_dis_max, P_max, P_ref, eta_ch, eta_dis, xi);
disp("Iteration time: "+toc+" s");
x = sol.x; % State of charge for each PEV at each time step
P = sol.P_ch-sol.P_dis; % Power for each PEV at each time step

%% Plots
% State of charge
figure, hold on, grid on;
for p = randi(m, 1, 10) % Plot random 10 PEVs
    plot(0:T:N*T, x(p, :), 'LineWidth', 1.5); % Plot state of charge over time for each PEV
end
title('\textbf{PEV state of charge}', 'Interpreter', 'LaTeX'); 
xlabel('Time [$h$]', 'Interpreter', 'LaTeX'); 
ylabel('Charge [$kWh$]', 'Interpreter', 'LaTeX'); 
xlim([0 N*T]), ylim([0 13]);
 
% Aggregated power
P_agg = sum(P, 1); % Calculate aggregated power at each time step
figure, hold on, grid on;
stairs(0:T:(N-1)*T, P_agg, 'LineWidth', 1.5); % Plot aggregated power over time
stairs(0:T:(N-1)*T, P_max, '--', 'LineWidth', 1.5); % Plot maximum power over time
stairs(0:T:(N-1)*T, P_ref, ':', 'LineWidth', 2); % Plot reference power over time
title('\textbf{Aggregated power profile}', 'Interpreter', 'LaTeX'); 
xlabel('Time [$h$]', 'Interpreter', 'LaTeX'); 
ylabel('Power [$kW$]', 'Interpreter', 'LaTeX'); 
xlim([0 (N-1)*T-0.01]), ylim([0 520]);
legend('Aggregated power', 'Maximum power', 'Reference power', 'Location', 'south', 'Interpreter', 'LaTeX');