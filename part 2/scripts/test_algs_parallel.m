clear, clc, close all;
tol = 1e-3;
rng(202405);

%% Parameters
N = 24;
T = 1/3;
m = 50; % multiple of 10
k_max = 20;

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

%% Parameters
P_tilde_star = NaN(N, m);
t = NaN(N, k_max+1);
rho = NaN(N, k_max);

%% Algorithm 1
k = 1;
rho(:, k) = zeros(N, 1);
J_tilde_star = inf;

while(k <= k_max)
    fprintf('%d) ', k)
    l_max = 300;
    alpha = 0.01./((1:l_max));
    P = NaN(N, m, l_max+1);
    P_hat = NaN(N, m, l_max+1);
    P_LP_hat = NaN(N, m);
    P_lambda_star = NaN(N, m);

    lambda = NaN(N, l_max+1);
    h = 2;
    l = 2;
    lambda(:, l) = zeros(N, 1);
    count = 0;
    while(l <= l_max)
        parfor p = 1:m
            pevs(p) = pevMpcIter(pevs(p), x_init(p), F(p), lambda(:, l), zeros(N, 1), zeros(N, 1));
            P(:, p, l) = pevs(p).sol.P;
        end
        lambda(:, l+1) = max(zeros(N, 1), lambda(:, l) + alpha(l)*(sum(P(:, :, l), 2) - P_max + rho(:, k)));

        if l >= 2 && max(abs(lambda(:, l+1) - lambda(:, l))) < 0.01 && max(abs(lambda(:, l) - lambda(:, l-1))) < 0.01
            count = count + 1;
            fprintf('|');
            parfor p = 1:m
                den = sum(alpha((l-h+1):l));
                num = sum(reshape(alpha((l-h+1):l), 1, 1, []).*P(:, p, (l-h+1):l), 3);
                P_hat(:, p, h) = num./den;
            end

            if h >= 3 && max(abs(P_hat(:, :, h)-P_hat(:, :, h-1)), [], 'all') < 0.1 && max(abs(P_hat(:, :, h-1)-P_hat(:, :, h-2)), [], 'all') < 0.1
                fprintf('!');
                P_LP_hat = P_hat(:, :, h);

                P_target = P(:, :, (l-h+1):l);
                parfor p = 1:m
                    to_minimize = NaN(size(P_target, 3), 1);
                    for j = 1:size(to_minimize, 1)
                        to_minimize(j) = dot(pevs(p).xi, P_target(:, p, j));
                    end
                    [~, ix] = min(to_minimize);
                    P_lambda_star(:, p) = P_target(:, p, ix);
                end
                break;
            end
            h = h + 1;
        end
        l = l + 1;
    end

    J = sum(dot([pevs.xi], P_lambda_star, 1), 2); %+ sum(t(:, k), 1);
    fprintf(' J = %f', J);

    P_agg = sum(P_lambda_star, 2);

    figure(1), hold off;
    stairs(0:T:(N-1)*T, P_agg, 'LineWidth', 1.5);
    hold on, grid on;
    stairs(0:T:(N-1)*T, P_max, '--', 'LineWidth', 1.5);
    stairs(0:T:(N-1)*T, rho(:, k), 'r', 'LineWidth', 1.5);
    title('\textbf{Aggregated power profile}', 'Interpreter', 'LaTeX');
    xlabel('Time [$h$]', 'Interpreter', 'LaTeX');
    ylabel('Power [$kW$]', 'Interpreter', 'LaTeX');
    xlim([0 (N-1)*T-0.01]), ylim([0 max(P_max)+10]);
    legend('Aggregated power', 'Maximum power', '$\rho$', 'Location', 'north', 'Interpreter', 'LaTeX');
    drawnow;

    margin = min(ones(N, 1)-P_agg./P_max)*100;
    fprintf(' - margin: %.1f %%', margin);
    if all(P_agg <= P_max + tol) && (J < J_tilde_star)
        P_tilde_star = P_lambda_star;
        J_tilde_star = J;
        fprintf(' <-- best feasible solution so far');
        beep;
    end

    fprintf('\n');
    rho(:, k+1) = max(zeros(N, 1), sum(P_lambda_star-P_LP_hat, 2));
    % rho(:, k+1) = sum(P_lambda_star-P_LP_hat, 2);
    % rho(:, k+1) = max(zeros(N, 1), rho(:, k) + beta(k)*(sum(P_lambda_star-P_LP_hat, 2) - rho(:, k)));
    % rho(:, k+1) = max(rho(:, k), rho(:, k) + sum(P_lambda_star-P_LP_hat, 2));
    % rho(:, k+1) = max(zeros(N, 1), rho(:, k) + sum(P_lambda_star-P_LP_hat, 2));
    % rho(:, k+1) = 100./P_max;
    % rho(:, k+1) = rho(:, k) + sum(P_lambda_star, 2) - P_max;
    % rho(:, k+1) = rho(:, k) + sum(P_lambda_star-P_LP_hat, 2);
    % rho2(:, k+1) = max(zeros(N, 1), rho2(:, k) + sum(P_lambda_star, 2) - P_ref);
    % rho3(:, k+1) = max(zeros(N, 1), rho3(:, k) + sum(P_lambda_star, 2) + P_ref);

    k = k + 1;
end