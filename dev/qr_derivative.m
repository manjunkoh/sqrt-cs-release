% Script to reproduce the results from "Differentiating the QR Decomposition"

clear; clc;

% Set the random number generator seed for reproducibility.
% Note: MATLAB and R use different random number generators, so the numerical
% results will not be identical to the paper's, but the quality of the
% approximation will be analogous.
rng(12345, 'twister');

% Create random matrices X and Y (the perturbation dX), as in the example.
% X is 10x3 and Y is a small perturbation.
x = randn(10, 3);
y = randn(10, 3) / 100;

% Call the function to compute the QR decompositions and their differentials.
[qx, rx, qz, rz, dq, dr] = d_qr(x, y);


%% --- Results Comparison ---

% Error of the zero-order approximation (i.e., just using Qx and Rx).
% This corresponds to "Qz - Qx" and "Rz - Rx".
err_Q0 = sum(abs(qz - qx), 'all');
err_R0 = sum(abs(rz - rx), 'all');

% Error of the first-order (linear) approximation.
% This corresponds to "Qz - (Qx + dQ)" and "Rz - (Rx + dR)".
err_Q1 = sum(abs(qz - (qx + dq)), 'all');
err_R1 = sum(abs(rz - (rx + dr)), 'all');

% --- Display the results in a format similar to the paper ---
fprintf('--- Reproduction of Paper''s Example ---\n\n');
fprintf('Paper''s Results:\n');
fprintf('  Sum of absolute diffs for Qz - Qx:         0.1316906\n');
fprintf('  Sum of absolute diffs for Rz - Rx:         0.0698008\n');
fprintf('  Sum of absolute diffs for Qz - (Qx + dQ):  0.0014707\n');
fprintf('  Sum of absolute diffs for Rz - (Rx + dR):  0.0014424\n\n');

fprintf('MATLAB Results:\n');
fprintf('  Sum of absolute diffs for Qz - Qx:         %.7f\n', err_Q0);
fprintf('  Sum of absolute diffs for Rz - Rx:         %.7f\n', err_R0);
fprintf('  Sum of absolute diffs for Qz - (Qx + dQ):  %.7f\n', err_Q1);
fprintf('  Sum of absolute diffs for Rz - (Rx + dR):  %.7f\n', err_R1);

fprintf('\nConclusion:\n');
fprintf('The error for the linear approximation (Qx+dQ, Rx+dR) is about\n');
fprintf('two orders of magnitude smaller than the zero-order approximation,\n');
fprintf('demonstrating the quality of the derived differentials.\n');