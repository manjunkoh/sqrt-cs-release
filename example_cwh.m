%% Example: CWH (Clohessy-Wiltshire-Hill) covariance steering
% Relative orbital motion control using CWH equations

clear; clc;
addpath ./SCvxStar/src/
addpath(genpath('./utils'))
addpath('./src')

figure_settings

flag_solved_qr = false;

%% System Parameters 
% Chief orbit parameters
r0 = 7228;  % Chief orbit radius (km) - altitude ~850 km
mu_earth = 3.986004418e5;  % Earth's gravitational parameter (km^3/s^2)

% Time discretization
dt = 30;  % Time step (sec)
N = 14;   % Number of intervals
tof = dt*N;
t_his = linspace(0, tof, N+1);

% Create CWH dynamical system
DS = astro.CWH(r0, mu_earth);

% System dimensions
nx = DS.nx;  % 6 (position + velocity in 3D)
nu = DS.nu;  % 3 (acceleration control in 3D)
nw = nx;      % Process noise dimension after discretization

% Uncertainty Parameters
% Initial dispersion
sigma_r0_m = 100;    % Initial position dispersion (m)
sigma_v0_mps = 1.0;    % Initial velocity dispersion (m/s)

% Terminal dispersion
sigma_rf_m = 10.0;   % Terminal position std dev (m)
sigma_vf_mps = 0.1;    % Terminal velocity std dev (m/s)

% Stochastic acceleration
sigma_a_mmps32 = 1.0;     % Stochastic acceleration (mm/s^(3/2))

% Convert units: initial parameters are in m, but we work in km
% State is in km and km/s, so we need to convert
sigma_r0_km = sigma_r0_m / 1000;  % 0.1 km
sigma_v0_kms = sigma_v0_mps / 1000; % 0.001 km/s
sigma_rf_km = sigma_rf_m / 1000;  % 0.01 km
sigma_vf_kms = sigma_vf_mps / 1000; % 0.0001 km/s
sigma_a_kms32 = sigma_a_mmps32 / 1E6;  % 1e-6 km/s^(3/2)

% Initial covariance
P0 = zeros(nx, nx);
P0(1:3, 1:3) = sigma_r0_km^2 * eye(3);  % Position covariance
P0(4:6, 4:6) = sigma_v0_kms^2 * eye(3); % Velocity covariance

% Terminal covariance
P_f = zeros(nx, nx);
P_f(1:3, 1:3) = sigma_rf_km^2 * eye(3);  % Position covariance
P_f(4:6, 4:6) = sigma_vf_kms^2 * eye(3); % Velocity covariance

% Mean Boundary Conditions 
% Initial mean state: r_0 = [-3.0, 0.126, 0] km, v_0 = [0, 0, 0] km/s
mu_0 = [-3.0; 0.126; 0; 0; 0; 0];  % [r; v] in km and km/s

% Terminal mean state: r_f = [0.0, 0.05, 0]^T km, v_f = [0, 0, 0] km/s
mu_f = [0.0; 0.05; 0; 0; 0; 0];    % [r; v] in km and km/s

g = @(t,x) sigma_a_kms32 * [zeros(3); eye(3)];
[A_sys, B_sys, ~, G_sys] = DS.discretize_LT_SDE(zeros(nx,N), zeros(nu,N), t_his, g);

% Comment out for low thrust
% for k = 1:N
    % B_sys(:,:,k) = A_sys(:,:,k) * [zeros(3); eye(3)];
% end
% B_sys = repmat([zeros(3); eye(3)], [1, 1, N]);
% A_sys = cat(3, A_sys, eye(6));
% B_sys = cat(3, B_sys, [zeros(3); eye(3)]);
% G_sys = cat(3, G_sys, zeros(nx, nx));
% N = N + 1;

% Chance Constraints
% Risk bounds
epsilon_x = 1e-3;  % State risk bound
epsilon_u = 1e-3;  % Control risk bound

% Control constraint: max ΔV magnitude = 10.0 m/s
u_max_ms = 10.0;  % m/s
u_max_kms = u_max_ms / 1000;  % Convert to km/s

% Control chance constraint: P(||u||_2 <= u_max) >= 1 - epsilon_u
control_chance_constraint = {
    struct('type', 'norm', 'gamma', u_max_kms, 'p', epsilon_u, 'n', nu, 'nodes', 1:N)
};

% Note: State chance constraints not specified, but could be added
state_chance_constraint = {};

% Objective Weights for LQG objective
Q = 0.01 * eye(nx);  % State cost (small since we care mainly about terminal covariance)
R = 1 * eye(nu);   % Control cost

%% Solve with FullCovarianceSteering
disp('=== Solving CWH Problem with FullCovarianceSteering ===');

% Reference covariance for linearization of objective and chance constraints
% P_ref = 1.5 * eye(nx) * max([sigma_r0_km^2, sigma_v0_kms^2]);
Y_ref = 2E-6 * eye(nu);

% Optionally try scaling the covariance dynamics.
d = 1;

prob_fc = FullCovarianceSteering(...
    A=A_sys, B=B_sys, G=G_sys, ...
    P_0=P0, P_f=P_f, ...
    Q=Q, R=R, ...
    N=N, ...
    objective_type='LQG',...
    chance_constraints_control=control_chance_constraint, ...
    Y_ref = Y_ref, ...
    covariance_scaling = d, ...
    mu_0=mu_0, mu_f=mu_f);

tic

sdp_settings = sdpsettings('verbose', 2, 'solver', 'mosek', 'savesolveroutput',1, 'savesolverinput', 1);

% dump MOSEK data as a text file 
% sdp_settings.mosektaskfile = 'data/mosek_dump.ptf';


% sdp_settings.mosek.MSK_IPAR_LOG_INTPNT = 10;

% See link for tips on debugging numerical issues: https://docs.mosek.com/11.0/toolbox/debugging-numerical.html

% 1. MOSEK automatically chooses whether to solve the primal or dual.
% Manually set this by choosing 'MSK_SOLVE_PRIMAL' or 'MSK_SOLVE_DUAL'
% sdp_settings.mosek.MSK_IPAR_INTPNT_SOLVE_FORM = 'MSK_SOLVE_DUAL';

% 2. MOSEK automatically chooses whether to presolve the problem. 
% Manually set this by choosing 'MSK_PRESOLVE_MODE_ON' or 'MSK_PRESOLVE_MODE_OFF'
% sdp_settings.mosek.MSK_IPAR_PRESOLVE_USE = 'MSK_PRESOLVE_MODE_ON';

% 3. MOSEK by default uses the maximum number of available threads
% Manually set the number of threads.
% sdp_settings.mosek.MSK_IPAR_NUM_THREADS = 1;

diagnostic_fc = prob_fc.solve(sdp_settings);
time_fc = toc;

if diagnostic_fc.problem == 0
    fprintf('FullCovarianceSteering solved successfully in %.3f seconds\n', time_fc);
    fprintf('  Objective value: %.6f\n', prob_fc.optimal_objective);
    prob_fc.check_lossless();
    prob_fc.control_covariance_is_SDP(true);

elseif diagnostic_fc.problem == 4
    fprintf('FullCovarianceSteering experienced numerical issues\n')
    fprintf('  MOSEK error message: %s\n', diagnostic_fc.solveroutput.res.rcodestr)
    fprintf('  Objective value: %.6f\n', prob_fc.optimal_objective);
    prob_fc.check_lossless();
    prob_fc.control_covariance_is_SDP(true);

else
    disp('FullCovarianceSteering solver failed.');
    disp(yalmiperror(diagnostic_fc.problem));
    prob_fc = [];
end

%% Solve with SqrtQRCovarianceSteering
disp('=== Solving CWH Problem with SqrtQRCovarianceSteering ===');

% Initial guess - Linear interpolation between initial and final states
S0 = chol(P0, 'lower');
S_f = chol(P_f, 'lower');
init_guess = struct();
init_guess.S = interpolate_lower_triangular(S0, S_f, N+1, 'log-cholesky');
init_guess.L = zeros(nu, nx, N);
init_guess.mu = linspace_vec(mu_0, mu_f, N+1);
init_guess.v = zeros(nu, N);

% K_init = - dlqr(A_sys(:,:,1), B_sys(:,:,1), Q, R);
% 
% for k = 1:N
%     init_guess.L(:,:,k) = K_init * init_guess.S(:,:,k);
% end

% Forward propagation of initial covariance
% P_k = P0;
% init_guess.S(:,:,1) = chol(P_k, 'lower');
% for k = 1:N
%     A = A_sys(:,:,k);
%     B = B_sys(:,:,k);
%     G = G_sys(:,:,k);
% 
%     P_k = (A+B*K_init) * P_k * (A+B*K_init)' + G * G';
%     init_guess.S(:,:,k+1) = chol(P_k, 'lower');
%     init_guess.L(:,:,k) = K_init * init_guess.S(:,:,k);
% end

% P_k = P_f;
% init_guess.S(:,:,N+1) = chol(P_k, 'lower');
% for k = (N+1):-1:floor(N/2)+1 
%     A = A_sys(:,:,k-1);
%     B = B_sys(:,:,k-1);
%     G = G_sys(:,:,k-1);
% 
%     P_k = inv(A+B*K_init) * (P_k - G * G') * inv(A+B*K_init)';
%     init_guess.S(:,:,k-1) = chol(P_k, 'lower');
%     init_guess.L(:,:,k-1) = K_init * init_guess.S(:,:,k-1);
% end

% for k = 1:N
%     init_guess.L(:,:,k) = K_init * init_guess.S(:,:,k);
% end
% figure
% plot_traj_with_cov_ellipses(init_guess.mu, init_guess.S, init_guess=true);


prob_qr = SqrtQRCovarianceSteering(init_guess, ...
    N=N, ...
    A_sys=A_sys, B_sys=B_sys, G_sys=G_sys, ...
    P_0=P0, P_f=P_f, Q=Q, R=R, ...
    objective_type='DV99', ...
    mu_0=mu_0, mu_f=mu_f);

prob_qr.D = diag([10, 10, 10, 1000, 1000, 1000]);

scp_params = SCPParams();
scp_params.tol_opt = 1E-5;
scp_params.tol_feas = 1E-5;
scp_params.w_init = 0.1;
% scp_params.alpha1 = 1.5;
% scp_params.alpha2 = 2;
% scp_params.penalty_method = 'AL_p-norm';
% scp_params.w_init = 10;
% scp_params.r_min = 1E-10;
% scp_params.r_init = 1.0;

flag_solved_qr = prob_qr.solve(save_bool=false, scp_params=scp_params);
time_qr = seconds(prob_qr.scp.report.time);
J_qr = NaN;

if flag_solved_qr
    prob_qr.postprocess();
    % Compute objective value (LQG objective)
    J_qr = 0;
    for k = 1:prob_qr.N
        J_qr = J_qr + trace(Q * prob_qr.P(:,:,k)) + trace(R * prob_qr.P_u(:,:,k)) ...
            + prob_qr.mu(:,k)' * Q * prob_qr.mu(:,k) + prob_qr.v(:,k)' * R * prob_qr.v(:,k);
    end
    fprintf('SqrtQRCovarianceSteering solved successfully in %.3f seconds\n', time_qr);
    fprintf('  Number of iterations: %d\n', prob_qr.scp.report.iters);
    fprintf('  Objective value: %.6f\n', J_qr);
else
    disp('SqrtQRCovarianceSteering solver failed. Displaying information from last iteration...');
    disp(prob_qr.scp.this_iter)
    prob_qr.scp.plot_iter_history()
end

%% Monte Carlo Simulation for QR Method
if flag_solved_qr && ~isempty(prob_qr)
    fprintf('\n=== Running Monte Carlo Simulation for SqrtQRCovarianceSteering ===\n');
    
    % Number of Monte Carlo samples
    num_samples = 100;
    
    % Initialize storage for Monte Carlo trajectories
    x_mc = zeros(nx, N+1, num_samples);  % State trajectories
    u_mc = zeros(nu, N, num_samples);    % Control trajectories
    
    % Generate random initial states from N(mu_0, P0)
    S0_chol = chol(P0, 'lower');
    rng(42);  % Set seed for reproducibility
    x0_samples = mu_0 + S0_chol * randn(nx, num_samples);
    
    % Simulate forward for each sample
    for i = 1:num_samples
        x_k = x0_samples(:, i);
        x_mc(:, 1, i) = x_k;
        
        for k = 1:N
            % Compute control: u_k = v_k + K_k * (x_k - mu_k)
            u_k = prob_qr.v(:, k) + prob_qr.K(:,:,k) * (x_k - prob_qr.mu(:, k));
            u_mc(:, k, i) = u_k;
            
            % Generate process noise: w_k ~ N(0, I)
            w_k = randn(nx, 1);
            
            % Forward dynamics: x_{k+1} = A_k * x_k + B_k * u_k + G_k * w_k
            A_k = A_sys(:,:,k);
            B_k = B_sys(:,:,k);
            G_k = G_sys(:,:,k);
            x_k = A_k * x_k + B_k * u_k + G_k * w_k;
            x_mc(:, k+1, i) = x_k;
        end
    end
    
    % Compute Monte Carlo statistics
    mu_mc = mean(x_mc, 3);  % Mean trajectory
    P_mc = zeros(nx, nx, N+1);
    for k = 1:N+1
        P_mc(:,:,k) = cov(squeeze(x_mc(:,k,:))');
    end
    
    % Check terminal constraint satisfaction
    x_final = squeeze(x_mc(:, N+1, :));  % Final states (nx x num_samples)
    % Check if terminal mean is within tolerance
    mu_final_mc = mean(x_final, 2);
    mu_error = norm(mu_final_mc - mu_f);
    fprintf('  Terminal mean error: %.6f km (target: %.6f km)\n', mu_error, norm(mu_f));
    
    % Check terminal covariance constraint: P_final <= P_f (in the sense of eigenvalues)
    P_final_mc = cov(x_final');
    % Check if largest eigenvalue of P_final is less than largest eigenvalue of P_f
    lambda_max_final = max(eig(P_final_mc));
    lambda_max_target = max(eig(P_f));
    fprintf('  Terminal covariance: max eigenvalue = %.6e (target: %.6e)\n', lambda_max_final, lambda_max_target);
    
    % Check Mahalanobis distance: fraction of samples within 3-sigma of target
    P_f_inv = inv(P_f);
    mahal_dist = zeros(num_samples, 1);
    for i = 1:num_samples
        dx = x_final(:, i) - mu_f;
        mahal_dist(i) = sqrt(dx' * P_f_inv * dx);
    end
    within_3sigma = sum(mahal_dist <= 3) / num_samples;
    fprintf('  Samples within 3-sigma of target: %.1f%%\n', within_3sigma * 100);
    
    % Check control constraint violations (if applicable)
    if ~isempty(control_chance_constraint)
        u_norm_samples = vecnorm(u_mc, 2, 1);  % (N x num_samples)
        u_norm_max = max(u_norm_samples(:)) * 1000;  % Convert to m/s
        u_violations = sum(u_norm_samples(:) > u_max_kms) / (N * num_samples) * 100;
        fprintf('  Control constraint violations: %.2f%% (max ||u|| = %.3f m/s, limit = %.1f m/s)\n', ...
            u_violations, u_norm_max, u_max_ms);
    end
    
    fprintf('Monte Carlo simulation completed.\n');
end

%% Plot Results
% Define flags for which solutions to plot
use_fc = exist('prob_fc', 'var') && ~isempty(prob_fc) && (diagnostic_fc.problem == 0 || diagnostic_fc.problem == 4) && prob_fc.check_lossless();
use_qr = exist('prob_qr', 'var') && flag_solved_qr && ~isempty(prob_qr);

%%
if use_fc || use_qr
    
    figure(Position=[0, 0, 20, 10])
    hold on
    
    % Plot FullCovarianceSteering results
    if use_fc
        % Plot covariance ellipses in XY plane
        for k = 1:N+1
            P_pos_fc = prob_fc.P([1,2], [1,2], k);
            plot3sigmaEllipse(prob_fc.mu([1,2], k), P_pos_fc, 'g', 'HandleVisibility', 'off');
        end
        % Plot mean trajectory (XY plane)
        plot(prob_fc.mu(1,:), prob_fc.mu(2,:), 'g+-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'FullCov');
        % Plot mean control vectors
        quiver2d(prob_fc.mu(1:2,1:N), prob_fc.v(1:2,1:N), 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Nominal control');
    end
    
    % Plot SqrtQRCovarianceSteering results
    if use_qr
        % Plot covariance ellipses in XY plane (every 3rd step)
        for k = 1:N+1
            P_pos_qr = prob_qr.P([1,2], [1,2], k);
            plot3sigmaEllipse(prob_qr.mu([1,2], k), P_pos_qr, 'Color', [.6, .6, .6], 'HandleVisibility', 'off');
        end

        % Plot Monte Carlo trajectories (if available)
        if exist('x_mc', 'var')
            % Plot a subset of Monte Carlo trajectories (for clarity)
            num_samples_mc = size(x_mc, 3);
            num_traj_plot = min(100, num_samples_mc);
            plot_indices = round(linspace(1, num_samples_mc, num_traj_plot));
            
            for i = plot_indices
                plot(x_mc(1,:,i), x_mc(2,:,i), '.-', 'MarkerSize', 8, 'Color', [0.7, 0.85, 1.0, 0.5], 'LineWidth', 0.5, 'DisplayName', 'MC samples');
            end
            
            % Plot Monte Carlo mean trajectory
            % if exist('mu_mc', 'var')
            %     plot(mu_mc(1,:), mu_mc(2,:), 'b--', 'LineWidth', 1.2, 'DisplayName', 'Monte Carlo mean');
            % end
            
            % Plot final Monte Carlo samples
            % if exist('x_final', 'var')
            %     scatter(x_final(1,:), x_final(2,:), 15, 'b', 'filled', 'MarkerFaceAlpha', 0.2, 'DisplayName', 'Final samples');
            % end
        end

        % Plot mean trajectory (XY plane)
        plot(prob_qr.mu(1,:), prob_qr.mu(2,:), 'Color', [.6, .6, .6], 'Marker', '+', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '($\mu_k, P_k$)');
        % Plot mean control vectors
        quiver2d(prob_qr.mu(1:2,1:N), prob_qr.v(1:2,1:N), 'Color', '#0082B2', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Nominal control $v_k$');
        
    end
    
    % Plot initial and terminal conditions
    % plot3sigmaEllipse(mu_0([1,2]), P0([1,2], [1,2]), 'b', 'LineWidth', 2, 'DisplayName', 'Initial ($\mu_{\mathrm{init}}, P_{\mathrm{init}}$)');
    plot3sigmaEllipse(mu_f([1,2]), P_f([1,2], [1,2]), '--', 'Color', '#D55E00', 'LineWidth', 2, 'DisplayName', 'Terminal ($\mu_{\mathrm{fin}}, P_{\mathrm{fin}}$)');
    
    xlabel('$x$ (km)', 'Interpreter', 'latex')
    ylabel('$y$ (km)', 'Interpreter', 'latex')

    if exist('x_mc', 'var') && use_qr
        legend(legendUnq(), 'Location', 'south', 'NumColumns', 2, 'EdgeColor', 'white', 'IconColumnWidth', 20)
    else
        legend(legendUnq(), 'Location', 'south', 'NumColumns', 1, 'EdgeColor', 'white', 'IconColumnWidth', 20)
    end
    grid on
    axis equal

    add_zoomed_axis(gca, [-0.05, 0.05, 0, 0.10], [-0.6, 0.4, 0.52, 0.52])
    
    % exportgraphics(gcf, 'figures/example_cwh_trajectory.png', Resolution=300)
    % exportgraphics(gcf, 'figures/example_cwh_trajectory.pdf', ContentType='vector')
    
    %% Plot Control History
    figure(Position=[0, 0, 15, 17])
    tiledlayout(3, 1)
    
    % Plot u_x, u_y components
    component_labels = {'$u_x$', '$u_y$'};
    for comp_idx = 1:2
        nexttile
        hold on
        
        % Collect all y-values for this component to set appropriate y-limits
        y_values = [];
        
        if use_fc
            v_comp = prob_fc.v(comp_idx,:) * 1000;  % Convert to m/s
            % Compute 3-sigma bounds
            sigma_comp = 3 * sqrt(squeeze(prob_fc.P_u(comp_idx,comp_idx,:))) * 1000;  % Convert to m/s
            upper_comp = v_comp + sigma_comp';
            lower_comp = v_comp - sigma_comp';
            
            y_values = [y_values, v_comp, upper_comp, lower_comp];
            
            stairsZOH(t_his, upper_comp, 'g:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            stairsZOH(t_his, lower_comp, 'g:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            stairsZOH(t_his, v_comp, 'g-', 'LineWidth', 1.5, 'DisplayName', 'FullCov');
        end
        
        if use_qr
            v_comp_qr = prob_qr.v(comp_idx,:) * 1000;  % Convert to m/s
            % Compute 3-sigma bounds
            sigma_comp_qr = 3 * sqrt(squeeze(prob_qr.P_u(comp_idx,comp_idx,:))) * 1000;  % Convert to m/s
            upper_comp_qr = v_comp_qr + sigma_comp_qr';
            lower_comp_qr = v_comp_qr - sigma_comp_qr';
            
            y_values = [y_values, v_comp_qr, upper_comp_qr, lower_comp_qr];
            
            % Get stair-step coordinates for upper and lower bounds (without plotting yet)
            [t_plot, upper_comp_qr_plot] = stairs(t_his, [upper_comp_qr, upper_comp_qr(end)]);
            [~,  lower_comp_qr_plot] = stairs(t_his, [lower_comp_qr, lower_comp_qr(end)]);
            % Create filled region between upper and lower bounds
            % For stairs plot: go forward along upper, then backward along lower
            patch_x = [t_plot(:); flipud(t_plot(:))];
            patch_y = [upper_comp_qr_plot(:); flipud(lower_comp_qr_plot(:))];
            % Plot patch first (so it's behind the lines)
            patch(patch_x, patch_y, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '$3\sigma$ bounds');

            % Add Monte Carlo control data if available
            if exist('u_mc', 'var')
                % Plot a subset of individual Monte Carlo control trajectories
                num_samples_mc = size(u_mc, 3);
                num_traj_plot = min(100, num_samples_mc);
                plot_indices = round(linspace(1, num_samples_mc, num_traj_plot));
                
                u_mc_comp = squeeze(u_mc(comp_idx, :, :)) * 1000;  % Convert to m/s (N x num_samples)
                
                % Collect MC y-values
                y_values = [y_values, u_mc_comp(:)'];
                
                % Plot individual trajectories first (so they appear behind other elements)
                for i = plot_indices
                    stairsZOH(t_his, u_mc_comp(:, i)', 'Color', [0.7, 0.85, 1.0], 'LineWidth', 0.5, 'DisplayName', 'MC samples', 'LineStyle', '-');
                end
                
                % Compute mean and std of Monte Carlo controls for this component
                % u_mc_mean = mean(u_mc_comp, 2)';  % Mean across samples (1 x N)
                % u_mc_std = std(u_mc_comp, 0, 2)';  % Std across samples (1 x N)
                % u_mc_upper = u_mc_mean + 3 * u_mc_std;
                % u_mc_lower = u_mc_mean - 3 * u_mc_std;
                % 
                % Plot Monte Carlo mean and bounds
                % [t_plot_mc, upper_mc_plot] = stairs(t_his, [u_mc_upper, u_mc_upper(end)]);
                % [~, lower_mc_plot] = stairs(t_his, [u_mc_lower, u_mc_lower(end)]);
                % patch_x_mc = [t_plot_mc(:); flipud(t_plot_mc(:))];
                % patch_y_mc = [upper_mc_plot(:); flipud(lower_mc_plot(:))];
                % patch(patch_x_mc, patch_y_mc, [0.7, 0.85, 1.0], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'MC $3\sigma$ bounds');
                % stairsZOH(t_his, u_mc_upper, 'Color', [0.7, 0.85, 1.0], 'LineStyle', ':', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                % stairsZOH(t_his, u_mc_lower, 'Color', [0.7, 0.85, 1.0], 'LineStyle', ':', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                % stairsZOH(t_his, u_mc_mean, 'Color', [0.7, 0.85, 1.0], 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'MC mean');
            end

            % Now plot the bounds and mean
            % stairsZOH(t_his, upper_comp_qr, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            % stairsZOH(t_his, lower_comp_qr, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            stairsZOH(t_his, v_comp_qr, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Nominal control $v_k$');

        end
        
        % Set y-limits with margin
        if ~isempty(y_values)
            y_min = min(y_values);
            y_max = max(y_values);
            y_range = y_max - y_min;
            margin = 0.1 * y_range;  % 10% margin
            if y_range == 0
                margin = 0.1;  % Default margin if no range
            end
            ylim([y_min - margin, y_max + margin]);
        end
        
        ylabel([component_labels{comp_idx}, ' (m/s$^2$)'], 'Interpreter', 'latex')
        grid on
        xlim([t_his(1), t_his(end)])
        % if comp_idx == 2
            legend(legendUnq(), 'Location', 'north', 'EdgeColor', 'white', 'IconColumnWidth', 20)
        % end

    end
    
    % Plot control norm
    nexttile
    hold on
    
    % Collect all y-values for control norm to set appropriate y-limits
    y_values_norm = [];
    
    if use_fc
        u_norm_fc = vecnorm(prob_fc.v, 2, 1) * 1000;  % Convert to m/s
        % For norm, approximate 3-sigma using largest eigenvalue of P_u
        sigma_norm_fc = zeros(1, N);
        for k = 1:N
            lambda_max_k = lambda_max(prob_fc.P_u(:,:,k));
            sigma_norm_fc(k) = sqrt(chi2inv(0.99, nu)) * lambda_max_k * 1000;  % Convert to m/s
        end
        upper_norm_fc = u_norm_fc + sigma_norm_fc;
        lower_norm_fc = max(0, u_norm_fc - sigma_norm_fc);  % Norm can't be negative
        
        y_values_norm = [y_values_norm, u_norm_fc, upper_norm_fc, lower_norm_fc];
        
        % Get stair-step coordinates for upper and lower bounds
        [t_plot_fc, upper_norm_fc_plot] = stairs(t_his, [upper_norm_fc, upper_norm_fc(end)]);
        [~, lower_norm_fc_plot] = stairs(t_his, [lower_norm_fc, lower_norm_fc(end)]);
        % Create filled region between upper and lower bounds
        patch_x_fc = [t_plot_fc(:); flipud(t_plot_fc(:))];
        patch_y_fc = [upper_norm_fc_plot(:); flipud(lower_norm_fc_plot(:))];
        % Plot patch first (so it's behind the lines)
        patch(patch_x_fc, patch_y_fc, 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '$3\sigma$ bounds');
        % Now plot the bounds and mean
        stairsZOH(t_his, upper_norm_fc, 'g:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        stairsZOH(t_his, lower_norm_fc, 'g:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        stairsZOH(t_his, u_norm_fc, 'g-', 'LineWidth', 1.5, 'DisplayName', 'FullCov');
        % Add constraint line
        % yline(u_max_ms, 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('Max: %.1f m/s', u_max_ms));
    end
    if use_qr
        u_norm_qr = vecnorm(prob_qr.v, 2, 1) * 1000;  % Convert to m/s
        % For norm, approximate 3-sigma using largest eigenvalue of P_u
        sigma_norm_qr = zeros(1, N);
        for k = 1:N
            lambda_max_k = lambda_max(prob_qr.P_u(:,:,k));
            sigma_norm_qr(k) = sqrt(chi2inv(0.99, nu)) * sqrt(lambda_max_k) * 1000;  % Convert to m/s
        end
        upper_norm_qr = u_norm_qr + sigma_norm_qr;
        lower_norm_qr = max(0, u_norm_qr - sigma_norm_qr);  % Norm can't be negative
        
        y_values_norm = [y_values_norm, u_norm_qr, upper_norm_qr, lower_norm_qr];
        
        % Get stair-step coordinates for upper and lower bounds
        [t_plot_qr, upper_norm_qr_plot] = stairs(t_his, [upper_norm_qr, upper_norm_qr(end)]);
        [~, lower_norm_qr_plot] = stairs(t_his, [lower_norm_qr, lower_norm_qr(end)]);
        % Create filled region between upper and lower bounds
        patch_x_qr = [t_plot_qr(:); flipud(t_plot_qr(:))];
        patch_y_qr = [upper_norm_qr_plot(:); flipud(lower_norm_qr_plot(:))];
        % Plot patch first (so it's behind the lines)
        patch(patch_x_qr, patch_y_qr, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '$99\%$ bounds');

        % Add Monte Carlo control norm if available
        if exist('u_mc', 'var')
            % Compute norm for each Monte Carlo sample
            u_norm_mc = zeros(N, size(u_mc, 3));
            for i = 1:size(u_mc, 3)
                u_norm_mc(:, i) = vecnorm(u_mc(:, :, i), 2, 1)' * 1000;  % Convert to m/s
            end
            
            % Collect MC y-values
            y_values_norm = [y_values_norm, u_norm_mc(:)'];
            
            % Plot a subset of individual Monte Carlo control norm trajectories
            num_samples_mc = size(u_mc, 3);
            num_traj_plot = min(100, num_samples_mc);
            plot_indices = round(linspace(1, num_samples_mc, num_traj_plot));
            
            % Plot individual trajectories first (so they appear behind other elements)
            for i = plot_indices
                stairsZOH(t_his, u_norm_mc(:, i)', 'Color', [0.7, 0.85, 1.0], 'LineWidth', 0.5, 'DisplayName', 'MC samples', 'LineStyle', '-');
            end
            
            u_norm_mc_mean = mean(u_norm_mc, 2)';  % Mean across samples (1 x N)
            % u_norm_mc_std = std(u_norm_mc, 0, 2)';  % Std across samples (1 x N)
            % u_norm_mc_upper = u_norm_mc_mean + 3 * u_norm_mc_std;
            % u_norm_mc_lower = max(0, u_norm_mc_mean - 3 * u_norm_mc_std);  % Norm can't be negative
            % 
            % % Plot Monte Carlo norm mean and bounds
            % [t_plot_mc, upper_norm_mc_plot] = stairs(t_his, [u_norm_mc_upper, u_norm_mc_upper(end)]);
            % [~, lower_norm_mc_plot] = stairs(t_his, [u_norm_mc_lower, u_norm_mc_lower(end)]);
            % patch_x_mc = [t_plot_mc(:); flipud(t_plot_mc(:))];
            % patch_y_mc = [upper_norm_mc_plot(:); flipud(lower_norm_mc_plot(:))];
            % patch(patch_x_mc, patch_y_mc, [0.7, 0.85, 1.0], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'MC $3\sigma$ bounds');
            % stairsZOH(t_his, u_norm_mc_upper, 'Color', [0.7, 0.85, 1.0], 'LineStyle', ':', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            % stairsZOH(t_his, u_norm_mc_lower, 'Color', [0.7, 0.85, 1.0], 'LineStyle', ':', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            % stairsZOH(t_his, u_norm_mc_mean, 'Color', [0.7, 0.85, 1.0], 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'MC mean');
        end

        % Now plot the bounds and mean
        % stairsZOH(t_his, upper_norm_qr, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        % stairsZOH(t_his, lower_norm_qr, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        stairsZOH(t_his, u_norm_qr, 'k--', 'LineWidth', 1.5, 'DisplayName', '$\|v_k\|$');

    end
    
    % Set y-limits with margin for control norm
    if ~isempty(y_values_norm)
        y_min_norm = min(y_values_norm);
        y_max_norm = max(y_values_norm);
        y_range_norm = y_max_norm - y_min_norm;
        margin_norm = 0.1 * y_range_norm;  % 10% margin
        if y_range_norm == 0
            margin_norm = 0.1;  % Default margin if no range
        end
        % Apply the same margin at both top and bottom
        y_min_lim = y_min_norm - margin_norm;  % Same margin at bottom
        y_max_lim = y_max_norm + margin_norm;  % Same margin at top
        ylim([y_min_lim, y_max_lim]);
    end
    
    xlabel('Time (s)', 'Interpreter', 'latex')
    ylabel('$\|u\|_2$ (m/s$^2$)', 'Interpreter', 'latex')
    % legend('Location', 'southeast')
    grid on
    xlim([t_his(1), t_his(end)])
    legend(legendUnq(), 'Location', 'north', 'EdgeColor', 'white', 'IconColumnWidth', 20)


    % exportgraphics(gcf, 'figures/example_cwh_control_history.png', Resolution=300)
    % exportgraphics(gcf, 'figures/example_cwh_control_history.pdf', ContentType='vector')
    
end

%% Summary
fprintf('\n=== Solution Summary for CWH Covariance Steering ===\n');
if (diagnostic_fc.problem == 0 || diagnostic_fc.problem == 4) && ~isempty(prob_fc)
    fprintf('FullCovarianceSteering: Time = %.3f s, Objective = %.6f\n', time_fc, prob_fc.optimal_objective);
else
    fprintf('FullCovarianceSteering: Failed\n');
end

if flag_solved_qr && ~isempty(prob_qr)
    fprintf('SqrtQRCovarianceSteering: Time = %.3f s, Objective = %.6f, Iterations = %d\n', time_qr, J_qr, prob_qr.scp.report.iters);
else
    fprintf('SqrtQRCovarianceSteering: Failed\n');
end

%% Compare Covariance Propagation Loss
fprintf('\n=== Covariance Propagation Loss Comparison ===\n');
if diagnostic_fc.problem == 0 || diagnostic_fc.problem == 4
    try
        loss_fc = prob_fc.compute_covariance_propagation_loss();
        fprintf('FullCovarianceSteering:\n');
        fprintf('  Mean loss: %.6e\n', mean(loss_fc));
        fprintf('  Max loss:  %.6e\n', max(loss_fc));
        fprintf('  Min loss:  %.6e\n', min(loss_fc));
    catch ME
        fprintf('FullCovarianceSteering: Error computing loss - %s\n', ME.message);
        loss_fc = [];
    end
else
    loss_fc = [];
end

if use_qr
    try
        loss_qr = prob_qr.compute_covariance_propagation_loss();
        fprintf('SqrtQRCovarianceSteering:\n');
        fprintf('  Mean loss: %.6e\n', mean(loss_qr));
        fprintf('  Max loss:  %.6e\n', max(loss_qr));
        fprintf('  Min loss:  %.6e\n', min(loss_qr));
    catch ME
        fprintf('SqrtQRCovarianceSteering: Error computing loss - %s\n', ME.message);
        loss_qr = [];
    end
else
    loss_qr = [];
end

% Plot comparison if both are availableCLose asdfasdfadf
if ~isempty(loss_fc) && ~isempty(loss_qr)
    figure(Position=[0, 0, 18, 8])
    hold on
    plot(0:N-1, loss_fc, 'o-', 'Color', '#0082B2', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Liu et al. (2025)');
    plot(0:N-1, loss_qr, 's-', 'Color', '#000000', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Proposed method');
    xlabel('Time step $k$', 'Interpreter', 'latex')
    ylabel('Loss')
    legend('Location', 'none', 'EdgeColor', 'white')
    yscale log

    ax = gca;
    ylims = ylim(ax);
    % Generate ticks at powers of 10 and intermediate values
    log_min = floor(log10(ylims(1)));
    log_max = ceil(log10(ylims(2)));
    yticks_val = 10.^(log_min:log_max);
    yticks(ax, yticks_val)

    exportgraphics(gcf, 'figures/example_cwh_covariance_propagation_loss.png', Resolution=300)
    exportgraphics(gcf, 'figures/example_cwh_covariance_propagation_loss.pdf', ContentType='vector')
end

%% Plot iteration history for QR method

prob_qr.scp.plot_iter_history(fig=figure(Position=[0, 0, 20, 15]), plot_delta=false)

% exportgraphics(gcf, 'figures/example_cwh_scp_iter_history.png', Resolution=300)
% exportgraphics(gcf, 'figures/example_cwh_scp_iter_history.pdf', ContentType='vector')