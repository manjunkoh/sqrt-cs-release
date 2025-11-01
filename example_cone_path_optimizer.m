%% Example: Cone-shaped path planning problem (from figure)
% Implements the system and constraints described in the attached figure.
% Uses BlockCholeskySteering (Okamoto et al. 2019 method)
clear; 
clc;
addpath ../SCvxStar/src/
addpath(genpath('./utils'))
addpath('./src')
yalmip('clear')

figure_settings

% Problem parameters
dt = 0.2;
N = 20;

% State: x = [x, vx, y, vy]'
nx = 4;
nu = 2;
nw = 4; % process noise dimension same as state here

% Build A and B (time-invariant)
A = [1  0   dt  0;
     0  1   0   dt;
     0  0   1   0;
     0  0   0   1];

B = [0.5*dt^2  0;
     0         0.5*dt^2;
     dt        0;
     0         dt];

% Replicate for time-varying arrays
A_sys = repmat(A, [1,1,N]);
B_sys = repmat(B, [1,1,N]);

% Process noise (diagonal)
G = diag([0.01, 0.01, 0.01, 0.01]);
G_sys = repmat(G, [1,1,N]);

% Covariance boundary conditions
Sigma0 = diag([0.1, 0.1, 0.01, 0.01]);
SigmaN = 0.5 * Sigma0;

% Mean boundary conditions
mu0 = [-10; 1; 0; 0];
muN = zeros(nx,1);

% Objective weights
Q = diag([10, 10, 1, 1]);
R = diag([1e3, 1e3]);

% Cone-shaped chance constraints: -0.2*(x-1) <= y <= 0.2*(x-1)
% Rearranged to alpha' * x <= beta form for state = [x, vx, y, vy]':
% Constraint 1: y <= 0.2(x-1)   ->  y <= 0.2x - 0.2  ->  0.2x + y <= 0.2
% Constraint 2: y >= -0.2(x-1)  ->  y >= -0.2x + 0.2  ->  0.2x - y <= 0.2
alpha1 = [0.2; 0; 1; 0]; beta1 = 0.2;  % Upper bound: y <= 0.2(x-1)
alpha2 = [0.2; 0; -1; 0]; beta2 = 0.2; % Lower bound: y >= -0.2(x-1)

% Probabilities for chance constraints
p_state = 0.005;

% Build chance constraint structs (apply to all time steps k=1:N+1)
state_cc = {
    struct('type', 'affine', 'alpha', alpha1, 'beta', beta1, 'p', p_state, 'nodes', 1:N+1);
    struct('type', 'affine', 'alpha', alpha2, 'beta', beta2, 'p', p_state, 'nodes', 1:N+1);
};

% Create problem using BlockCholeskySteering
prob = BlockCholeskySteering(...
    A=A_sys, B=B_sys, G=G_sys, ...
    P_0=Sigma0, P_f=SigmaN, ...
    Q=Q, R=R, ...
    N=N, ...
    chance_constraints_state=state_cc, ...
    mu_0=mu0, mu_f=muN);

% Solve
disp('Solving Block Cholesky Covariance Steering problem...')
tic
diagnostic = prob.solve('verbose', 0, 'solver', 'mosek');
toc

if diagnostic.problem == 0 || diagnostic.problem == 4
    % Get results
    P = prob.P;
    disp('Problem solved successfully!');
else
    disp('Solver failed.');
    disp(yalmiperror(diagnostic.problem));
    return
end

%% Plot results
figure
hold on

% Plot constraint boundaries (cone walls)
% State is [x, vx, y, vy]', so constraint alpha'*state <= beta becomes:
% alpha(1)*x + alpha(2)*vx + alpha(3)*y + alpha(4)*vy <= beta
% For plotting in (x,y) space, we set vx=0, vy=0 and solve for y:
% alpha(1)*x + alpha(3)*y = beta -> y = (beta - alpha(1)*x) / alpha(3) if alpha(3) != 0
x_ch = linspace(-10, 1, 100)';
for j = 1:length(state_cc)
    cc = state_cc{j};
    if strcmp(cc.type, 'affine')
        a = cc.alpha;
        b = cc.beta;
        if abs(a(3)) > 1e-10  % Check that constraint involves y (a(3) != 0)
            y_ch = (b - a(1)*x_ch) / a(3);
            plot(x_ch, y_ch, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Constraint boundary');
        end
    end
end

% Plot covariance ellipses along trajectory
for k = 1:N+1
    plot3sigmaEllipse(prob.mu(:,k), P(:,:,k), 'b')
end

% Plot initial and terminal constraints
plot3sigmaEllipse(mu0, Sigma0, 'r', 'LineWidth', 2, 'DisplayName', 'Initial');
plot3sigmaEllipse(muN, SigmaN, 'r--', 'LineWidth', 2, 'DisplayName', 'Terminal');

% Plot mean trajectory
plot(prob.mu(1,:), prob.mu(3,:), 'k+-', 'LineWidth', 1.5, 'DisplayName', 'Mean trajectory');

axis equal
xlabel("$x$", 'Interpreter', 'latex')
ylabel("$y$", 'Interpreter', 'latex')
title('Block Cholesky Covariance Steering with Cone-Shaped Chance Constraints')
legend('Location', 'best')
grid on
% exportgraphics(gcf, './figures/example_cone_path.png')