% Test propagation of square root covariance using QR decomposition

% clc; clear;

A = [1 0.2; 0 1];
D = [0.4 0; 0.4 0.6];
P_0 = [5 -1; -1 1];

nx = size(A,1);
nw = size(D,2);

S_0 = chol(P_0, 'lower');

S_1 = qr([A * S_0, D]', "econ")';

% if need the Q matrix as well, use:
[Q, S_1] = qr([A * S_0, D]', "econ");
S_1 = S_1';
(Q * S_1')' - [A * S_0, D]

assert( norm(S_1 * S_1' - (A * P_0 * A' + D * D')) < 1e-10 );