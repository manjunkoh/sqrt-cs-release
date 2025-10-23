%% Example: Quadrotor 2D path planning (parameters from the figure)
% Discrete triple integrator lateral/longitudinal dynamics
clear; clc;

% Time horizon and discretization
N = 60;
dt = 0.1;

% Dimensions
nx = 6; % triple integrator in x and y (position and velocity and accel?)
nu = 2; % two inputs (x and y thrust)
nw = 6;

% Build block matrices as in the figure
I2 = eye(2);
Z2 = zeros(2);
A_block = [I2, dt*I2, 0*I2;
           Z2, I2,  dt*I2;
           Z2, Z2,  I2];
% A_sys is time-varying; replicate
A_sys = repmat(A_block, [1,1,N]);

B_block = [Z2; Z2; dt*I2];
B_sys = repmat(B_block, [1,1,N]);

G_block = 0.1 * eye(nx); % small process noise
G_sys = repmat(G_block, [1,1,N]);

% Covariance boundary conditions
Sigma_i = eye(nx);
Sigma_f = 0.1 * eye(nx);

% Mean boundary conditions
mu_i = [20; 0; zeros(nx-2,1)];
mu_f = zeros(nx,1);

% Chance constraint parameters (boxes): example linear constraints
% p_x = 0.05;
% alpha_x = {[1;0;zeros(nx-2,1)], [-1;0;zeros(nx-2,1)], [0;1;zeros(nx-2,1)], [0;-1;zeros(nx-2,1)]};
% beta_x = [22, -3, 7, -7];

% p_u = 0.05;
% alpha_u = {[1;0], [-1;0], [0;1], [0;-1]};
% beta_u = [25, -25, 25, -25];

% Create affine chance constraints for state and control at specified times
% state_cc = cell(1, N+1);
% control_cc = cell(1, N);

% Example: apply alpha_x/beta_x at time steps 20 and 40
% for k = 1:N+1
%     state_cc{k} = {};
% end
% for k = [20, 40]
%     for j = 1:length(alpha_x)
%         a = alpha_x{j}; b = beta_x(j);
%         state_cc{k}{end+1} = AffineChanceConstraint(a, b, p_x);
%     end
% end
% 
% % Control constraints for all k
% for k = 1:N
%     control_cc{k} = {};
%     for j = 1:length(alpha_u)
%         a = alpha_u{j}; b = beta_u(j);
%         control_cc{k}{end+1} = AffineChanceConstraint(a, b, p_u);
%     end
% end

% Initial guess struct (S and L)
init.S = linspace_mat(chol(Sigma_i, 'lower'), chol(Sigma_f, 'lower'), N+1);
init.L = zeros(nu, nx, N);
init.mu = zeros(nx, N+1);
init.v = zeros(nu, N);

% Create problem
prob = SqrtQRCovarianceSteering(init, ...
    N=N, nx=nx, nu=nu, nw=nw, ...
    A_sys=A_sys, ...
    B_sys=B_sys, ...
    G_sys=G_sys, ...
    P_0=Sigma_i, P_f=Sigma_f, Q=eye(nx), R=eye(nu), ...
    objective_type='LQG', ...
    mu_0=mu_i, mu_f=mu_f);

flag_solved = prob.solve();

if flag_solved
    prob.postprocess();
    disp('Solved successfully.');
else
    disp('Solver failed.');
end
