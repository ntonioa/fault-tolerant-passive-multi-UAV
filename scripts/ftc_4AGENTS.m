clear, clc, close all;

%% graph
N = 4;
A = [0, 1, 0, 1;
    1, 0, 1, 0;
    0, 1, 0, 1;
    1, 0, 1, 0];
G = graph(A);
E = full(incidence(G));
L = full(laplacian(G));

param.N = N;
param.E = E;
param.A = A;
param.L = L;

%% parameters
M1 = 0.16*eye(2);
M = repmat(M1, 1, 1, N);
param.M = M;

B1 = 1.5*eye(2);
B = repmat(B1, 1, 1, N);
param.B = B;

K1 = 50*eye(2);
K = repmat(K1, 1, 1, N);
param.K = K;

R = eye(2);
param.R = R;

param.theta = 0.1;
param.tau_star = 2;

pos_init(:, 1) = [0; 0];
pos_init(:, 2) = [5; 0];
pos_init(:, 3) = [5; 5];
pos_init(:, 4) = [0; 5];
param.pos_init = pos_init;

d = zeros(2, N);
for j = 1:N
    for i = 1:N
        d(:, j) = d(:, j) + E(i, j)*pos_init(:, i);
    end
end
param.d = d;

csi_init = d;
param.csi_init = csi_init;

param.x_init = 15*ones(N, 1);

param.fault_i = 4; %2;
param.fault_perc = 0; %95; 
param.fault_time = 4;
param.dist_amp = 40; %0;

% figure;
% plot(G, 'XData', pos_init(1, :), 'YData', pos_init(2, :), 'NodeLabel', {'Agent 1', 'Agent 2', 'Agent 3', 'Agent 4'});

%% animations
video = 'traj.avi';
% video = 'ener.avi';
v = VideoWriter(video, 'Motion JPEG AVI');
v.FrameRate = 30;
open(v);

param.T = 0.02;

model = 'ftc_4AGENTS_sim';
out = sim(model);
scr_siz = get(0,'ScreenSize') ;

f1 = figure(1);
set(gcf, 'Color', 'w');
set(f1, 'Position', [100, 100, 750, 600]);
hold on;
grid on;
axis equal;
xlabel('X position [m]');
ylabel('Y position [m]');
title('Agent trajectories');

f2 = figure(2);
set(gcf, 'Color', 'w');
set(f2, 'Position', [100, 100, 750, 600]);
hold on;
grid on;
xlabel('Agent');
ylabel('Tank energy [J]');
title('Tank energy evolution');
ylim([0 max([out.H_tank.Data(:)])+10]);
xlim([0.5 4.5]);
ax = gca;
ax.XTick = unique( round(ax.XTick) );

numFrames = size(out.pos.Data, 3);

traj = gobjects(1, N);
ag = gobjects(1, N);
tanks = gobjects(1, N);
links = gobjects(1, N);
colors = {[0, 112, 192]/255, [38, 38, 38]/255, [38, 38, 38]/255, [38, 38, 38]/255};
for i = 1:N
    figure(f2);
    tanks(i) = bar(NaN, NaN, 'FaceColor', colors{i});

    figure(f1);
    traj(i) = plot(NaN, NaN, ':', 'LineWidth', 1.5, 'Color', colors{i});
    ag(i) = plot(NaN, NaN, 'o', 'LineWidth', 1.5, 'Color', colors{i}, 'DisplayName', 'Followers');
    links(i) = plot(NaN, NaN, '', 'LineWidth', 0.5, 'Color', [116, 116, 116]/255);
end
input = quiver(NaN, NaN, NaN, NaN, 'LineWidth', 1, 'Color', colors{1});
ref = plot(NaN, NaN, 'LineWidth', 1.5, 'Color', [146, 208, 80]/255);
ref_tip = plot(NaN, NaN, 'pentagram', 'LineWidth', 1.5, 'Color', [146, 208, 80]/255, 'DisplayName', 'Reference');

% set(ag(1), 'DisplayName', 'Leader');
% leg = legend([ag(1), ag(3), ref_tip],'location', 'best');

for t = 1:numFrames
    for i = 1:N
        set(traj(i), 'XData', squeeze(out.pos.Data(1, i, max([1 t-20]):t)), 'YData', squeeze(out.pos.Data(2, i, max([1 t-20]):t)));
        set(ag(i), 'XData', squeeze(out.pos.Data(1, i, t)), 'YData', squeeze(out.pos.Data(2, i, t)));
        set(tanks(i), 'XData', i, 'YData', squeeze(out.H_tank.Data(1, i, t)))
        set(links(i), 'XData', [squeeze(out.pos.Data(1, i, t)); squeeze(out.pos.Data(1, mod(i, N)+1, t))], 'YData', [squeeze(out.pos.Data(2, i, t)); squeeze(out.pos.Data(2, mod(i, N)+1, t))]);
    end
    set(ref, 'XData', squeeze(out.pos1_ref.Data(1, 1, 1:t)), 'YData', squeeze(out.pos1_ref.Data(2, 1, 1:t)));
    set(ref_tip, 'XData', squeeze(out.pos1_ref.Data(1, 1, t)), 'YData', squeeze(out.pos1_ref.Data(2, 1, t)));
    set(input, 'XData', squeeze(out.pos.Data(1, 1, t)), 'YData', squeeze(out.pos.Data(2, 1, t)), 'UData', squeeze(out.f_ex.Data(1, 1, t))/10, 'VData', squeeze(out.f_ex.Data(2, 1, t))/10, 'MaxHeadSize', 3/norm(out.f_ex.Data(:, 1, t)));

    if t*param.T >= param.fault_time
        set(ag(param.fault_i), 'Color', 'r');
        set(tanks(param.fault_i), 'FaceColor', 'r');
        if out.j0.Data(t) == 0 && mod(t, 20) < 10
            set(ag(param.fault_i), 'Color', [38, 38, 38]/255);
            set(tanks(param.fault_i), 'FaceColor', [38, 38, 38]/255);
        end
        if out.j0.Data(t) ~= 0
            temp = 1:N;
            temp(temp == out.j0.Data(t)) = [];
            for i = 1:N-1
                set(links(temp(i)), 'XData', [squeeze(out.pos.Data(1, temp(i), t)); squeeze(out.pos.Data(1, temp(mod(i, N-1)+1), t))], 'YData', [squeeze(out.pos.Data(2, temp(i), t)); squeeze(out.pos.Data(2, temp(mod(i, N-1)+1), t))]);
                set(links(out.j0.Data(t)), 'XData', [], 'YData', []);
            end
        end
    end

    drawnow;

    frame = getframe(f1);
    % frame = getframe(f2);
    writeVideo(v, frame);
end
close(v);

%%
f3 = figure(3);
set(gcf, 'Color', 'w');
set(f3, 'Position', [100, 100, 750, 600]);
hold on;
grid on;
xlabel('Time [s]');
ylabel('Energy [J]');
title('Energy curves');

plot(out.H.Time, squeeze(sum(out.H.Data(1, :, :), 2)), 'LineWidth', 1.5, 'DisplayName', 'Total');
plot(out.H_ex.Time, squeeze(out.H_ex.Data), 'LineWidth', 1.5, 'DisplayName', 'External');
xline(param.fault_time, '--r', 'LineWidth', 1.5, 'DisplayName', 'Time of fault');
xline(6.08, '--', 'LineWidth', 1.5,  'DisplayName', 'Time of reconfiguration');
leg = legend('location', 'best');

%% energy function for the spring
% Define the function
f = @(x,y) (1/2) * (x - 1).^2 + (1/2) * (y - 1).^2 + (1 ./ (x.^2 + y.^2));

% Define the range for x and y
x = linspace(-5, 5, 100);
y = linspace(-5, 5, 100);

% Create a grid of x and y values
[X, Y] = meshgrid(x, y);

% Evaluate the function at each grid point
Z = f(X, Y);

% Create a surface plot
% figure;
% surf(X, Y, Z);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title('Surface plot of z = (1/2)(x-1)^2 + (1/2)(y-1)^2 + 1/(x^2 + y^2)');
% colorbar;
% zlim([0 10]);
% clim([0, 10]);

% Create a contour plot
f4 = figure(4);
set(gcf, 'Color', 'w');
set(f4, 'Position', [100, 100, 750, 600]);
hold on;
grid on;
contour(X, Y, Z, linspace(0, 10, 20));
xlabel('\xi_1 / d_1');
ylabel('\xi_2 / d_2');
title('Example of virtual spring storage function');
colorbar;
xlim([-2 4]);
ylim([-2 4]);
zlim([0 8]);
clim([0, 8]);
%% test
% figure;
% hold on;
% grid on;

% plot(out.H_pot.Time, squeeze(sum(out.H_pot.Data(1, :, :), 2)) + squeeze(sum(out.H_tank.Data(1, :, :), 2)), 'LineWidth', 1.5);
% plot(out.H_pot.Time, squeeze(sum(out.H_kin.Data(1, :, :), 2)), 'LineWidth', 1.5);

% plot(out.H_pot.Time, squeeze(sum(out.H_tank.Data(1, :, :), 2)), 'LineWidth', 1.5);