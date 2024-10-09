% function[P_tilde_star] = algorithm1(P_max, P_ref, pevs, N, m, k_max, F, x_init)

%% Parameters
P_tilde_star = zeros(N, m);
P = zeros(N, m, k_max);
t = zeros(N, k_max);
rho = zeros(N, k_max);

%% Algorithm 1
k = 1;
J_tilde_star = inf;

while(k <= k_max)
    [P_LP_hat, P_lambda_star] = algorithm2(rho(:, k), P_max, P_ref, pevs, N, m, F, x_init);
    P(:, :, k) = P_lambda_star;
    t(:, k) = abs(sum(P(:, :, k), 2) - P_ref);

    J = sum(dot([pevs.xi], P(:, :, k), 1), 2) + sum(t(:, k), 1)
    
    if all(sum(P(:, :, k), 2) <= P_max) && (J <= J_tilde_star) % *
        P_tilde_star = P(:, :, k);
        J_tilde_star = J;
        disp('new solution found!')
    end
    rho(:, k+1) = max(zeros(N, 1), sum(P(:, :, k)-P_LP_hat, 2));
    k = k + 1;
end

% end