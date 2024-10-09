clear, clc, close all;

%% graph
N = 2;
A = ones(N) - eye(N);
G = graph(A);
E = full(incidence(G));
L = full(laplacian(G));

param.N = N;
param.E = E;

figure;
plot(G)

%% parameters
M1 = eye(3);
M = repmat(M1, 1, 1, N);
param.M = M;

B1 = 0.01*eye(3);
B = repmat(B1, 1, 1, N);
param.B = B;

K1 = 10*eye(3);
K = repmat(K1, 1, 1, nchoosek(N, 2));
param.K = K;

%% test
x = rand(1, N);
p = rand(3, N);

r = zeros(1, 3, N);
for i= 1:N
    r(:, :, i) = (B(:, :, i)*M(:, :, i)\p(:, i)).'./repmat(x(:, i), 1, 3);
end