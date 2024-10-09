% function [P_LP_hat, P_lambda_star] = algorithm2(rho, P_max, P_ref, pevs, N, m, F, x_init)
%%
l_max = 100;
P = zeros(N, m, l_max+1);
P_hat = zeros(N, m, l_max+1);
P_LP_hat = zeros(N, m);
P_lambda_star = zeros(N, m);

lambda = zeros(N, l_max+1);
mu = zeros(N, l_max+1);
nu = zeros(N, l_max+1);
alpha = 0.001./((1:l_max+1));

%%
h = 1;
l = 1;
count = 0;
while(l <= l_max)
    for p = 1:m
        pevs(p) = pevMpcIter(pevs(p), x_init(p), F(p), lambda(:, l) , mu(:, l), nu(:, l));
        P(:, p, l) = pevs(p).sol.P;
    end
    lambda(:, l+1) = max(zeros(N, 1), lambda(:, l) + alpha(l)*(sum(P(:, :, l), 2) - P_max + rho));
    mu(:, l+1) = max(zeros(N, 1), mu(:, l) + alpha(l)*(sum(P(:, :, l), 2) - P_ref));
    nu(:, l+1) = max(zeros(N, 1), nu(:, l) - alpha(l)*(sum(P(:, :, l), 2) - P_ref));
    theta = [lambda; mu; nu];

    if max(abs(theta(:, l+1) - theta(:, l))) < 0.1
        count = count + 1;
        fprintf('%d) theta is at convergence\n', count);
        for p = 1:m
            den = sum(alpha(l-h+1:l+1));
            num = sum(reshape(alpha(l-h+1:l+1), 1, 1, []).*P(:, p, l-h+1:l+1), 3);
            P_hat(:, p, h+1) = num./den;
        end

        if max(abs(P_hat(:, :, h+1)-P_hat(:, :, h)), [], 'all') < 0.1
            disp('P_hat is at convergence');
            P_LP_hat = P_hat(:, :, h+1);

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

% end