%% Aggregated power
figure, hold on, grid on;
stairs(0:T:(N-1)*T, sum(P_tilde_star, 2), 'LineWidth', 1.5); % Plot aggregated power
stairs(0:T:(N-1)*T, P_max, '--', 'LineWidth', 1.5); % Plot maximum power limit
title('\textbf{Aggregated power profile}', 'Interpreter', 'LaTeX'); 
xlabel('Time [$h$]', 'Interpreter', 'LaTeX'); 
ylabel('Power [$kW$]', 'Interpreter', 'LaTeX'); 
xlim([0 (N-1)*T-0.01]), ylim([0 max(P_max)+10]);
legend('Aggregated power', 'Maximum power', 'Location', 'north', 'Interpreter', 'LaTeX');

%% Global constraint violation
figure;
stairs(0:k-1, [margin(1:k-1), 0], 'LineWidth', 1.5); % Plot maximum power constraint violation
yline(tol*100, '--r', 'LineWidth', 1.5); % Plot threshold
title('\textbf{Maximum aggregated power constraint violation}', 'Interpreter', 'LaTeX'); 
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('Violation [$\%$]', 'Interpreter', 'LaTeX'); 
xlim([0 k-1.01]), ylim([min([margin, 0])-20 max(margin)+20]);
legend('Maximum constraint violation', 'Threshold', 'Location', 'north', 'Interpreter', 'LaTeX');
grid on;

%% State of charge
figure, hold on, grid on;
rand_p = randi(m, 1, 10); % Randomly select 10 PEVs
for p = rand_p
    plot(0:T:N*T, pevs(p).sol.x, 'LineWidth', 1.5); % Plot state of charge for each selected PEV
end
title('\textbf{PEV state of charge}', 'Interpreter', 'LaTeX'); 
xlabel('Time [$h$]', 'Interpreter', 'LaTeX'); 
ylabel('Charge [$kWh$]', 'Interpreter', 'LaTeX'); 
xlim([0 N*T]), ylim([0 13]);

%% Norm of rho
figure, hold on, grid on;
stairs(0:k-1, [vecnorm(rho(:, 1:k-1), 2, 1), 0], 'r', 'LineWidth', 1.5);
title('\textbf{Norm of the variable $\rho$}', 'Interpreter', 'LaTeX'); 
xlabel('Iteration', 'Interpreter', 'LaTeX'); 
ylabel('$||\rho||$', 'Interpreter', 'LaTeX'); 
xlim([0 k-1.01]);

%% Lambda
% figure;
% plot(0:k-1, lambda(1:N-1, :)', 'LineWidth', 1.5); % Plot multipliers lambda
% title('\textbf{Multiplier $\lambda$}', 'Interpreter', 'LaTeX'); 
% xlabel('Iteration', 'Interpreter', 'LaTeX'); 
% ylabel('$\lambda$', 'Interpreter', 'LaTeX'); 
% xlim([0 k-1]);
% grid on;