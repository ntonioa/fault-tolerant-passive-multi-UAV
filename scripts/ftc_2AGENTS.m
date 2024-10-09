clear, clc, close all;

%% graph
N = 2;
A = ones(N) - eye(N);
G = graph(A);
E = full(incidence(G));
L = full(laplacian(G));

param.N = N;
param.E = E;

%% parameters
M1 = eye(3);
M = repmat(M1, 1, 1, N);
param.M = M;

B1 = 2*eye(3);
B = repmat(B1, 1, 1, N);
param.B = B;

K = 100*eye(3);
param.K = K;
param.d = [5; 0; 0];

%% initial conditions
pos_init(:, 1) = [2; 0; 0];
pos_init(:, 2) = [-2; 0; 0];
csi_init = pos_init(:, 1) - pos_init(:, 2);

figure;
plot(G, 'XData', pos_init(1, 1:2), 'YData', [0, 0], 'NodeLabel', {'Agent 1', 'Agent 2'});