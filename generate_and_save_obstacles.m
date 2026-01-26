clc; clear;
addpath(genpath('./utils'))
addpath ./src
addpath ./SCvxStar/src/
addpath ./obstacle_path_planning
figure_settings

function [obstacle_centers, obstacle_radii] = generate_random_environment(num_obstacles, obstacle_radius, wall_y_pos, mu_0, mu_f, Sigma_0, Sigma_f)
    is_valid = false;
    while ~is_valid
        obstacle_centers = [rand(num_obstacles, 1) * 10, rand(num_obstacles, 1) * (wall_y_pos - (-3)) + (-3)];
        obstacle_radii = obstacle_radius * ones(1, num_obstacles);
        is_valid = is_valid_environment(obstacle_centers, obstacle_radii, mu_0, mu_f, Sigma_0, Sigma_f);
    end
end

%% Parameters
basic_parameters;

num_obstacles = 3;

% Number of cases to generate (will keep only those that solve)
num_cases_to_generate = 50;
max_attempts_per_case = 5; % Maximum attempts to find a solvable case

% Storage for valid cases
valid_cases = struct();
valid_cases.obstacle_centers = {};
valid_cases.obstacle_radii = {};
valid_cases.x_opt = {};
valid_cases.u_opt = {};
valid_cases.case_number = [];

case_counter = 0;
total_attempts = 0;

fprintf('Generating %d solvable obstacle cases...\n', num_cases_to_generate);
fprintf('==========================================\n');

rng(0); % Set seed for reproducibility

while case_counter < num_cases_to_generate && total_attempts < max_attempts_per_case * num_cases_to_generate
    total_attempts = total_attempts + 1;
    
    if mod(total_attempts, 10) == 0
        fprintf('Attempt %d: Found %d/%d valid cases...\n', total_attempts, case_counter, num_cases_to_generate);
    end
    
    % Generate random environment
    [obstacle_centers, obstacle_radii] = generate_random_environment(num_obstacles, obstacle_radius, wall_y_pos, mu_0, mu_f, Sigma_0, Sigma_f);
    
    % Try to solve with deterministic method
    penalty_scalar_obstacle = 100;
    max_iters = 50;
    flag_deterministic = false;
    
    yalmip('clear')
    x = sdpvar(nx, num_nodes+1, 'full');
    u = sdpvar(nu, num_nodes, 'full');
    lambda = sdpvar(1, num_nodes, 'full');
    
    x_ref = linspace_vec(mu_0, mu_f, num_nodes+1);
    
    for iter = 1:max_iters
        constraints = [];
        objective = 0;
        
        for k = 1:num_nodes
            constraints = [constraints, x(:,k+1) == A * x(:,k) + B * u(:,k)];
            constraints = [constraints, u(:,k) <= u_max];
            constraints = [constraints, u(:,k) >= -u_max];
        end
        
        constraints = [constraints, x(:,1) == mu_0, x(:,num_nodes+1) == mu_f];
        constraints = [constraints, x(2,:) <= wall_y_pos];
        
        for k = 1:num_nodes
            for i = 1:size(obstacle_centers, 1)
                a = - (x_ref(pos_idx,k) - obstacle_centers(i,:)');
                b = - 0.5 * norm(x_ref(pos_idx,k) - obstacle_centers(i,:)')^2  + 0.5 * obstacle_radii(i)^2 - a' * x_ref(pos_idx,k);
                constraints = [constraints
                    a' * x(pos_idx,k) + b <= lambda(k)
                ];
            end
        end
        
        constraints = [constraints, lambda >= 0];
        
        for k = 1:num_nodes
            objective = objective + x(:,k)' * Q * x(:,k) + u(:,k)' * R * u(:,k) + penalty_scalar_obstacle * lambda(k);
        end
        
        sol = optimize(constraints, objective, sdpsettings('verbose', 0));
        
        if sol.problem
            % Problem infeasible, break and try next case
            break;
        end
        
        x_opt = value(x);
        u_opt = value(u);
        
        if norm(x_opt - x_ref) < 1e-3 && is_collision_free(x_opt, obstacle_centers, obstacle_radii)
            flag_deterministic = true;
            break;
        end
        
        x_ref = x_opt;
    end
    
    % If solved successfully, save the case
    if flag_deterministic
        case_counter = case_counter + 1;
        valid_cases.obstacle_centers{case_counter} = obstacle_centers;
        valid_cases.obstacle_radii{case_counter} = obstacle_radii;
        valid_cases.x_opt{case_counter} = x_opt;
        valid_cases.u_opt{case_counter} = u_opt;
        valid_cases.case_number(case_counter) = case_counter;
        
        if mod(case_counter, 10) == 0
            fprintf('✓ Found valid case %d/%d\n', case_counter, num_cases_to_generate);
        end
    end
end

fprintf('\n==========================================\n');
fprintf('Generation complete!\n');
fprintf('Total attempts: %d\n', total_attempts);
fprintf('Valid cases found: %d\n', case_counter);
if total_attempts > 0
    fprintf('Success rate: %.2f%%\n', 100 * case_counter / total_attempts);
end

if case_counter == 0
    fprintf('WARNING: No valid cases found! Try increasing max_attempts_per_case or adjusting parameters.\n');
    return;
end

% Save to mat file
save_filename = './data/obstacle_cases.mat';
fprintf('\nSaving to %s...\n', save_filename);

% Create data directory if it doesn't exist
if ~exist('./data', 'dir')
    mkdir('./data');
end

save(save_filename, 'valid_cases', 'num_obstacles', 'obstacle_radius', 'num_cases_to_generate', 'total_attempts');

fprintf('Saved %d valid cases to %s\n', case_counter, save_filename);
fprintf('\nFile contains:\n');
fprintf('  - valid_cases.obstacle_centers: cell array of obstacle centers (one per case)\n');
fprintf('  - valid_cases.obstacle_radii: cell array of obstacle radii (one per case)\n');
fprintf('  - valid_cases.x_opt: cell array of optimal state trajectories\n');
fprintf('  - valid_cases.u_opt: cell array of optimal control sequences\n');
fprintf('  - valid_cases.case_number: array of case numbers\n');
