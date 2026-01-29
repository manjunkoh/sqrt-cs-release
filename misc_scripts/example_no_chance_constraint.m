% Very simple comparison of FullCovarianceSteering and SqrtQRCovarianceSteering
% using the unconstrained example from Liu 2025
clc; clear;
addpath ./SCvxStar/src/
addpath(genpath("./utils"))
addpath('./src')
yalmip("clear")
figure_settings()

%% Problem parameters (from Liu 2025)
A = [1 0.2; 0 1];
B = [0.02 0.2]';
G = [0.4 0; 0.4 0.6];
P_0 = [5 -1; -1 1];
P_f = [0.5 -0.4; -0.4 2];
mu_0 = [30; -5];
mu_f = [0; 0];
Q = 0.5 * eye(2);
R = 1;
nx = 2;
nu = 1;
nw = size(G, 2);
N = 29;

A_sys = repmat(A, [1, 1, N]);
B_sys = repmat(B, [1, 1, N]);
G_sys = repmat(G, [1, 1, N]);

%% Solve covariance steering using FullCovarianceSteering class

full_cs = FullCovarianceSteering(...
	A=A_sys, B=B_sys, G=G_sys, ...
	P_0=P_0, P_f=P_f, ...
	Q=Q, R=R, ...
	N=N, mu_0=mu_0, mu_f=mu_f);

diagnostic = full_cs.solve();

if diagnostic.problem == 0 || diagnostic.problem == 4
	P = full_cs.P;
	Y = full_cs.P_u;
	K = full_cs.K;
	
	% Compute covariance part of objective
	J_cov = 0;
    
	for k = 1:N
		J_cov = J_cov + trace(Q*P(:,:,k)) + trace(R * Y(:,:,k));
	end
else
	disp('FullCovarianceSteering solve not successful')
	disp(yalmiperror(diagnostic.problem))
	return
end

%% Perform interpolation of covariances to generate initial guess
% try changing between 'log-cholesky' and 'cholesky'
init_guess_struct.S = interpolate_lower_triangular(chol(P_0, 'lower'), chol(P_f, 'lower'), N+1, 'log-cholesky');
init_guess_struct.L = zeros(nu, nx, N);
init_guess_struct.mu = linspace_vec(mu_0, mu_f, N+1);
init_guess_struct.v = zeros(nu, N);

sqrt_cs = SqrtQRCovarianceSteeringOptimizer(...
    init_guess_struct,...
	N=N, ...
	A_sys=repmat(A, [1, 1, N]), ...
	B_sys=repmat(B, [1, 1, N]), ...
	G_sys=repmat(G, [1, 1, N]), ...
	P_0=P_0, P_f=P_f, Q=Q, R=R, ...
	mu_0=mu_0, mu_f=mu_f);

scp_params = SCPParams();
scp_params.tol_opt = 1E-4;
scp_params.tol_feas = 1E-4;

sqrt_cs.solve(scp_params = scp_params);

sqrt_cs.postprocess();

%% Plot the results
figure(Position=[0, 0, 30, 15])

tiledlayout(1, 2)

nexttile
hold on
for k = 1:N+1
	plot3sigmaEllipse(full_cs.mu(:,k), full_cs.P(:,:,k), 'b')
end
plot3sigmaEllipse(mu_0, P_0, 'r')
plot3sigmaEllipse(mu_f, P_f, 'r--')
axis equal
xlabel("$x_1$")
ylabel("$x_2$")

nexttile
hold on
for k = 1:N
	plot3sigmaEllipse(sqrt_cs.mu(:,k), sqrt_cs.P(:,:,k), 'b')
end
plot3sigmaEllipse(mu_0, P_0, 'r')
plot3sigmaEllipse(mu_f, P_f, 'r--')
axis equal
xlabel("$x_1$")
ylabel("$x_2$")

% exportgraphics(gcf, './figures/sqrt_covariance_steering_result.png');

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
% exportgraphics(gcf, './figures/sqrt_covariance_steering_convergence.png');