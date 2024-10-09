clear, clc, close all;

%% Parameters
% Simulation
N = 24; % Prediction horizon
T = 1/3; % Time step duration
m = 5; % Number of PEVs
k_max = 10; % Maximum number of iterations

% Global
P_max_1 = 5.0*ones(9, 1);
P_max_2 = 2.0*ones(6, 1);
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

%% Parameters
P_tilde_star = zeros(N, m);
P1 = zeros(N, m, k_max+1);
t = zeros(N, k_max+1);
rho = zeros(N, k_max);

%% Algorithm 1
k = 1;
J_tilde_star = inf;

while(k <= k_max)
    l_max = 100;
    alpha = 0.05./((1:l_max));
    beta = 0.005./((1:l_max));
    P = zeros(N, m, l_max+1);
    P_hat = zeros(N, m, l_max+1);
    P_LP_hat = zeros(N, m);
    P_lambda_star = zeros(N, m);

    lambda = zeros(N, l_max+1);
    mu = zeros(N, l_max+1);
    nu = zeros(N, l_max+1);
    h = 2;
    l = 2;
    count = 0;
    while(l <= l_max)
        for p = 1:m
            pevs(p) = pevMpcIter(pevs(p), x_init(p), F(p), lambda(:, l) , mu(:, l), nu(:, l));
            P(:, p, l) = pevs(p).sol.P;
        end
        lambda(:, l+1) = max(zeros(N, 1), lambda(:, l) + alpha(l)*(sum(P(:, :, l), 2) - P_max + rho(:, k)));
        % mu(:, l+1) = max(zeros(N, 1), mu(:, l) + alpha(l)*(sum(P(:, :, l), 2) - P_ref));
        % nu(:, l+1) = max(zeros(N, 1), nu(:, l) - alpha(l)*(sum(P(:, :, l), 2) - P_ref));
        theta = [lambda; mu; nu];

        if max(abs(theta(:, l) - theta(:, l-1))) < 0.02
            count = count + 1;
            fprintf('%d) theta is at convergence\n', count);
            for p = 1:m
                den = sum(beta((l-h+1):l));
                num = sum(reshape(beta((l-h+1):l), 1, 1, []).*P(:, p, (l-h+1):l), 3);
                P_hat(:, p, h) = num./den;
            end

            if max(abs(P_hat(:, :, h)-P_hat(:, :, h-1)), [], 'all') < 0.1
                disp('P_hat is at convergence');
                P_LP_hat = P_hat(:, :, h);

                for p = 1:m
                    P_target = P(:, p, (l-h+1):l); % *
                    to_minimize = zeros(size(P_target, 3), 1);
                    for j = 1:size(to_minimize, 1)
                        to_minimize(j) = dot(pevs(p).xi, P_target(:, 1, j));
                    end
                    [~, ix] = min(to_minimize);
                    P_lambda_star(:, p) = P_target(:, :, ix);
                end
                break;
            end
            h = h + 1;
        end
        l = l + 1;
    end


    t(:, k) = abs(sum(P_lambda_star, 2) - P_ref);

    J = sum(dot([pevs.xi], P_lambda_star, 1), 2) %+ sum(t(:, k), 1)

    P_agg = sum(P_lambda_star, 2); % Aggregated power profile
    figure(1), hold off, grid on;
    stairs(0:T:(N-1)*T, P_agg, 'LineWidth', 1.5); % Plot aggregated power
    hold on;
    stairs(0:T:(N-1)*T, P_max, '--', 'LineWidth', 1.5); % Plot maximum power limit
    stairs(0:T:(N-1)*T, rho(:, k), 'r', 'LineWidth', 1.5); % Plot rho
    title('\textbf{Aggregated power profile}', 'Interpreter', 'LaTeX');
    xlabel('Time [$h$]', 'Interpreter', 'LaTeX');
    ylabel('Power [$kW$]', 'Interpreter', 'LaTeX');
    xlim([0 (N-1)*T-0.01]), ylim([0 20]);
    legend('Aggregated power', 'Maximum power', '$\rho$', 'Location', 'north', 'Interpreter', 'LaTeX');
    drawnow;

    if all(sum(P_lambda_star, 2) <= P_max) && (J <= J_tilde_star) % *
        P_tilde_star = P_lambda_star;
        J_tilde_star = J;
        disp('new solution found!')
    end
    % rho(:, k+1) = max(zeros(N, 1), sum(P_lambda_star-P_LP_hat, 2));
    rho(:, k+1) = max(zeros(N, 1), rho(:, k) + sum(P_lambda_star-P_LP_hat, 2));
    
    % rho(:, k+1) = 100./P_max;
    k = k + 1;
end