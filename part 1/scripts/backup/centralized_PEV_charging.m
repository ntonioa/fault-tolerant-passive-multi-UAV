clear, clc, close all;

%% parameters
% simulation
N = 24;
m = 500;
T = 1/3;

% global
P_max_1 = 500*ones(1, 10);
P_max_2 = 200*ones(1, 4);
P_max = [P_max_1, P_max_2, P_max_1];
P_ref = [290*ones(1, N-1), 0];

% local
F = N*ones(m, 1);
P_ch_max = 5*ones(m, 1);
P_ch_min = 1.3*ones(m, 1);
P_dis_max = 0*ones(m, 1);
P_dis_min = 0*ones(m, 1);
x_max = 8*(ones(m, 1)+rand(m, 1));
x_min = 1*ones(m, 1);
x_init = (0.2*ones(m, 1)+0.3*rand(m, 1)).*x_max;
x_ref = (0.55*ones(m, 1)+0.25*rand(m, 1)).*x_max;
eta_ch = 0.925*ones(m, 1)+0.06*rand(m, 1);
eta_dis = 1*ones(m, 1);
xi = 0.3*rand(m, N);

%% solution
tic
sol = CenMpc(N, T, m, x_min, x_max, x_init, x_ref, F, P_ch_min, P_ch_max, P_dis_min, P_dis_max, P_max, P_ref, eta_ch, eta_dis, xi);
disp("Iteration time: "+toc+" s");
x = sol.x;
P = sol.P_ch-sol.P_dis;

%% plots
% state
figure, hold on, grid on;
for p = randi(m, 1, 10)
    plot(0:T:N*T, x(p, :), 'LineWidth', 1.5);
end
title('PEV state charge'), xlabel('Time [h]'), ylabel('Charge [kWh]'); 
xlim([0 N*T]), ylim([0 13]);
 
% aggregated power
P_agg = sum(P, 1);
figure, hold on, grid on;
stairs(0:T:(N-1)*T, P_agg, 'LineWidth', 1.5);
stairs(0:T:(N-1)*T, P_max, '--', 'LineWidth', 1.5);
stairs(0:T:(N-1)*T, P_ref, ':', 'LineWidth', 2);
title('Power profiles'), xlabel('Time [h]'), ylabel('Power [kW]');
xlim([0 (N-1)*T]), ylim([0 520]);
legend('Aggregated power', 'Maximum power', 'Reference power', 'Location', 'south');