% Square root covariance steering with QR decomposition-based covariance propagation
clc; clear;
addpath ../SCvxStar/src/
addpath(genpath("../utils"))
yalmip("clear")
figure_settings()

%%
A = [1 0.2; 0 1];
B = [0.02 0.2]';
D = [0.4 0; 0.4 0.6];
P_0 = [5 -1; -1 1];
P_f = [0.5 -0.4; -0.4 2];
mu_0 = [30; -5];
mu_f = [0; 0];
Q = 0.5 * eye(2);
R = 1;
nx = 2;
nu = 1;
nw = size(D, 2);
N = 29;

%% Lossless relaxation: minimize the sum of trace
P = sdpvar(nx,nx,N+1);
U = sdpvar(nu,nx,N);
Y = sdpvar(nu,nu,N);
mu = sdpvar(nx,N+1);
v = sdpvar(nu, N);

J = 0;
J_cov = 0;
for k = 1:N
	J_cov = J_cov + trace(Q*P(:,:,k)) + trace(R * Y(:,:,k));
	J = J ...
		+ trace(Q*P(:,:,k)) ...
		+ trace(R * Y(:,:,k)) ...
		+ mu(:,k)'*Q*mu(:,k) ...
		+ v(:,k)'*R*v(:,k);
end

constraints = [P(:,:,1) == P_0
	P_f - P(:,:,end) >= 0
	mu(:,1) == mu_0
	mu(:,end) == mu_f];

for k = 1:N
	constraints = [constraints
		Y(:,:,k) >= 0
		P(:,:,k+1) == A*P(:,:,k)*A.' + A*U(:,:,k).'*B.' + B*U(:,:,k)*A.' + B*Y(:,:,k)*B.' + D * D.'
		[P(:,:,k) U(:,:,k).'; U(:,:,k) Y(:,:,k)] >= 0
		mu(:,k+1) == A*mu(:,k) + B*v(:,k)];
end

for k= 1:N+1
	constraints = [constraints
		P(:,:,k) >= 0];
end

options = sdpsettings('verbose', 0, 'solver', 'mosek');
diagnostic = optimize(constraints, J, options);

if diagnostic.problem == 0
	P = value(P);
	U = value(U);
	Y = value(Y);
	mu = value(mu);
	v = value(v);
    J = value(J);
    J_cov = value(J_cov);
else
	disp('YALMIP solve not successful')
	disp(yalmiperror(diagnostic.problem))
	return
end

% Retrieve feedback gains
K = zeros(nu, nx, N);
for i = 1:N
	K(:,:,i) = U(:,:,i) / P(:,:,i);
end

% Plot
figure
hold on
for k = 1:N+1
	plot3sigmaEllipse(mu(:,k), P(:,:,k), 'b')
end
plot3sigmaEllipse(mu_0, P_0, 'r')
plot3sigmaEllipse(mu_f, P_f, 'r--')
axis equal
xlabel("$x_1$")
ylabel("$x_2$")
exportgraphics(gcf, './figures/full_covariance_steering_result.png');

disp("Objective value: " + J)
disp("Covariance part of objective: " + J_cov)

%% Use the true solution as the initial guess -- just for testing
% init_guess_struct.S = zeros(nx, nx, N+1);
% init_guess_struct.L = zeros(nu, nx, N);
% for k = 1:N+1
% 	init_guess_struct.S(:,:,k) = chol(P(:,:,k), 'lower');
% end
% 
% for k = 1:N
% 	init_guess_struct.L(:,:,k) = K(:,:,k) * chol(P(:,:,k), 'lower');
% end

%% Perform interpolation of covariances to generate initial guess
S_0 = chol(P_0, 'lower');
S_f = chol(P_f, 'lower');
init_guess_struct.S = linspace_mat(S_0, S_f, N+1);
% this might be a bad init guess
init_guess_struct.L = ones(nu, nx, N);
init_guess_struct.mu = zeros(nx, N+1);
init_guess_struct.v = zeros(nu, N);
%%
sqrt_cs = SqrtQRCovarianceSteering(init_guess_struct,...
	N=N, nx=nx, nu=nu, nw=nw, ...
	A_sys=repmat(A, [1, 1, N]), ...
	B_sys=repmat(B, [1, 1, N]), ...
	G_sys=repmat(D, [1, 1, N]), ...
	P_0=P_0, P_f=P_f, Q=Q, R=R, ...
	mu_0=mu_0, mu_f=mu_f);

scp_params = SCPParams();
scp_params.tol_opt = 1E-4;
scp_params.tol_feas = 1E-6;

sqrt_cs.solve(scp_params = scp_params);

sqrt_cs.postprocess();

%% Plot the results
figure
hold on
for k = 1:N
	plot3sigmaEllipse(sqrt_cs.mu(:,k), sqrt_cs.P(:,:,k), 'b')
end
plot3sigmaEllipse(mu_0, P_0, 'r')
plot3sigmaEllipse(mu_f, P_f, 'r--')
axis equal
xlabel("$x_1$")
ylabel("$x_2$")
exportgraphics(gcf, './figures/sqrt_covariance_steering_result.png');
%%
J_cov_sqrt = 0;
for k = 1:N
	J_cov_sqrt = J_cov_sqrt + trace(Q*sqrt_cs.P(:,:,k)) ...
		+ trace(R * sqrt_cs.P_u(:,:,k));
end

disp("Covariance part of objective (sqrt method): " + value(J_cov_sqrt));
disp("Covariance part of objective (full covariance method): " + value(J_cov))

%%
sqrt_cs.scp.plot_iter_history()
exportgraphics(gcf, './figures/sqrt_covariance_steering_convergence.png');
