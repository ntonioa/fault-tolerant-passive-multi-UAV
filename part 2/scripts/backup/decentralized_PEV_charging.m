clear, clc, close all;

%% parameters
% simulation
N = 24;
m = 500;
T = 1/3;
k_max = 100;

% global
P_max_1 = 500*ones(1, 10);
P_max_2 = 200*ones(1, 4);
P_max = [P_max_1, P_max_2, P_max_1];
P_ref = [290*ones(1, N-1), 0];
tol = 30;

% local
F = N*ones(m, 1);
x_init = zeros(m, 1);
pevs(1:m, 1) = PevMpc;
for p = 1:m
    x_max = 8*(1+rand);
    x_init(p) = (0.2+0.3*rand)*x_max;
    x_ref = (0.55+0.25*rand)*x_max;
    eta_ch = 0.925+0.06*rand;
    xi = 0.3*rand(1, N);

    pevs(p, 1) = PevMpc(N, T, x_max, 1, x_ref, 5, 1.3, 0, 0, eta_ch, 1, xi);
end

%% variables
P = zeros(m, N, k_max+1);

rho = zeros(m, N, k_max);
rho_agg = zeros(N, k_max);
lambda = zeros(N, k_max);
mu = zeros(N, k_max);
nu = zeros(N, k_max);
alpha = 0.001./((2:k_max));

%% solution
figure, grid on;
title('Current power profiles'), xlabel('Time [h]'), ylabel('Power [kW]');
k = 0;
while(true)
    k = k+1;
    disp(" ");
    disp("Iteration "+(k-1));
    
    % parallel computing
    P_next = zeros(m, N);
    lambda_curr = lambda(:, k);
    mu_curr = mu(:, k);
    ni_curr = nu(:, k);
    time = 0;
    parfor p = 1:m
        tic;
        pevs(p) = pevMpcIter(pevs(p), x_init(p), F(p), lambda_curr, mu_curr, ni_curr);
        P_next(p, :) = pevs(p).sol.P;
        time = time+toc;
    end
    P(:, :, k+1) = P_next;
    disp("Average iteration time (per PEV): "+(time/m)+" s");
    % ------------------

    P_max_viol = max(sum(P(:, :, k+1), 1)-P_max);
    P_ref_gap = max(abs(sum(P(:, :, k+1), 1)-min(P_ref, P_max)));
    disp("Maximum violation of the maximum aggregated power constraint: "+P_max_viol+" kW");
    disp("Maximum deviation from the adjusted aggregated power reference: "+P_ref_gap+" kW");
    if  (P_max_viol <= 0 && P_ref_gap <= tol) || k == k_max 
        P = P(:, :, 1:k+1);
        rho = rho(:, :, 1:k);
        rho_agg = rho_agg(:, 1:k);
        lambda = lambda(:, 1:k);
        mu = mu(:, 1:k);
        nu = nu(:, 1:k);
        break;
    end

    rho(:, :, k+1) = reshape([pevs.s_up]-[pevs.s_down], m, []);
    rho_agg(:, k+1) = N*max(rho(:, :, k+1), [], 1);
    lambda(:, k+1) = max(zeros(N, 1), lambda(:, k)+alpha(k)*(sum(P(:, :, k+1), 1)'-P_max'+rho_agg(:, k+1)));
    mu(:, k+1) = max(zeros(N, 1), mu(:, k)+alpha(k)*(sum(P(:, :, k+1), 1)'-P_ref'));
    nu(:, k+1) = max(zeros(N, 1), nu(:, k)-alpha(k)*(sum(P(:, :, k+1), 1)'-P_ref'));

    % aggregated power step-by-step
    stairs(0:T:(N-1)*T, P_max, 'r');
    hold on;
    stairs(0:T:(N-1)*T, P_ref, '--');
    stairs(0:T:(N-1)*T, sum(P(:, :, k+1), 1));
    drawnow;
    hold off;
end

%% plots
% aggregated power
P_agg = sum(P(:, :, k+1), 1);
figure, hold on, grid on;
stairs(0:T:(N-1)*T, P_agg, 'LineWidth', 1.5);
stairs(0:T:(N-1)*T, P_max, '--', 'LineWidth', 1.5);
stairs(0:T:(N-1)*T, P_ref, ':', 'LineWidth', 2);
title('Power profiles'), xlabel('Time [h]'), ylabel('Power [kW]');
xlim([0 (N-1)*T]), ylim([0 520]);
legend('Aggregated power', 'Maximum power', 'Reference power', 'Location', 'south');

% global constraint violation
P_max_viol = squeeze(max(sum(P(:, :, 2:k+1), 1)-P_max));
P_ref_gap = squeeze(max(abs(sum(P(:, :, 2:k+1), 1)-min(P_ref, P_max))));
figure;

subplot(2, 1, 1);
plot(0:k-1, P_max_viol, 'LineWidth', 1.5);
yline(0, '--r', 'LineWidth', 1.5);
title('Maximum power global constraint violation'); 
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('Violation [$kW$]', 'Interpreter', 'LaTeX'); 
xlim([0 k-1]), ylim([min(P_max_viol)-20 max(P_max_viol)+20]);
legend('Maximum power violation', 'Threshold', 'Location', 'northeast');
grid on;

subplot(2, 1, 2);
plot(0:k-1, P_ref_gap, 'LineWidth', 1.5);
yline(tol, '--r', 'LineWidth', 1.5);
title('Maximum reference power gap'); 
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('Gap [$kW$]', 'Interpreter', 'LaTeX'); 
xlim([0 k-1]), ylim([min(P_ref_gap)-20 max(P_ref_gap)+20]);
legend('Maximum reference gap', 'Threshold', 'Location', 'northeast');
grid on;

% state
figure, hold on, grid on;
rand_p = randi(m, 1, 10);
for p = rand_p
    plot(0:T:N*T, pevs(p).sol.x, 'LineWidth', 1.5);
end
title('PEV state of charge'), xlabel('Time [h]'), ylabel('Charge [kWh]'); 
xlim([0 N*T]), ylim([0 13]);

% norm of rho
figure, hold on, grid on;
for p = rand_p
    plot(0:k-1, vecnorm(squeeze(rho(p, :, :))), 'LineWidth', 1.5);
end
title('Norm of the variable \rho'), xlabel('Time [h]'), ylabel('Charge [kWh]');
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('$||\rho||$', 'Interpreter', 'LaTeX'); 

% lambda
figure;
plot(0:k-1, lambda(1:N-1, :)', 'LineWidth', 1.5);
title('Multiplier \lambda'); 
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('$\lambda$', 'Interpreter', 'LaTeX'); 
grid on;

% mu
figure;
plot(0:k-1, mu(1:N-1, :)', 'LineWidth', 1.5);
title('Multiplier \mu'); 
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('$\mu$', 'Interpreter', 'LaTeX'); 
grid on;

% ni
figure;
plot(0:k-1, nu(1:N-1, :)', 'LineWidth', 1.5);
title('Multiplier \nu'); 
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('$\nu$', 'Interpreter', 'LaTeX'); 
grid on;
