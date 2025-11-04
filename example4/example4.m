%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 4: Lorenz system with control input to dot{y} term
% This example demonstrates trajectory optimization for the Lorenz system
% with a control input added to the y-dot equation.
% The implementation of the problem is in ExampleClass4.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')) % Path to SCvx* source code
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils')) % Path to utility functions

prob = ExampleClass4();
SCvxParams = SCPParams();
SCvxParams.tol_opt = 1E-3;
% SCvxParams.penalty_method = "AL_p-norm";
% SCvxParams.penalty_power = 1.1;
prob.solve(scp_params = SCvxParams, save_bool = false, clear_every_iter=true);

%% Propagate the trajectory using the converged control history
u_handle = @(t) getZOH(t, prob.t_his, prob.sol.u);
traj = prob.DS.propagate_with_LT(prob.x_init, [prob.t_his(1), prob.t_his(end)], u_handle);

%% Plot the trajectory in 3D
figure;
view(3)
hold on;
plot3d(prob.init_guess_struct.x, 'r--', DisplayName='Initial guess', LineWidth=1.5)
plot3d(traj.y(1:3,:), 'k', DisplayName='Optimal trajectory', LineWidth=2)

% Mark initial and final states
plot3(prob.x_init(1), prob.x_init(2), prob.x_init(3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Initial state')
plot3(prob.x_fin(1), prob.x_fin(2), prob.x_fin(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Final state')

xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('Lorenz System Trajectory')
legend()
axis equal
grid on
view(45, 20)
%% Plot individual state components vs time
figure;
subplot(2,1,1)
hold on
plot(prob.t_his, prob.sol.x(1,:), 'k-', DisplayName='$x$', LineWidth=2)
plot(prob.t_his, prob.sol.x(2,:), 'b-', DisplayName='$y$', LineWidth=2)
plot(prob.t_his, prob.sol.x(3,:), 'g-', DisplayName='$z$', LineWidth=2)
xlabel('time [s]')
ylabel('State')
title('State Trajectories')
legend(Orientation='horizontal')
grid on

subplot(2,1,2)
hold on
stairsZOH(prob.t_his, prob.sol.u, 'k', DisplayName='Control input $u$', LineWidth=1.5)
yline([prob.u_min, prob.u_max], 'r--', HandleVisibility='off')
xlabel('time [s]')
ylabel('Control input $u$')
title('Control Profile')
legend()
grid on

%% Plot the SCvx* convergence history
prob.scp.plot_iter_history()

%% Plot the processing time of the algorithm
prob.scp.plot_time()