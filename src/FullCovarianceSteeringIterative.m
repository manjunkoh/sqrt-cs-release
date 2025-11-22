classdef FullCovarianceSteeringIterative < handle
	% FullCovarianceSteeringIterative - Wrapper for FullCovarianceSteering
	% 
	% This class performs iterative solving of FullCovarianceSteering problems
	% with chance constraints by updating the reference solution (P_ref, Y_ref)
	% used for linearization.
	%
	% The algorithm:
	% 1. First solves the problem without chance constraints to obtain initial
	%    reference values (P_ref, Y_ref)
	% 2. Iteratively solves with chance constraints, updating P_ref and Y_ref
	%    from the previous solution until convergence
	%
	% Usage:
	%   wrapper = FullCovarianceSteeringIterative(...
	%       'A', A, 'B', B, 'G', G, 'P_0', P_0, 'P_f', P_f, ...
	%       'Q', Q, 'R', R, 'N', N, ...
	%       'chance_constraints_state', state_cc, ...
	%       'mu_0', mu_0, 'mu_f', mu_f);
	%   [diagnostic, prob] = wrapper.solve(sdp_settings);
	
	properties
		% Problem parameters (passed to FullCovarianceSteering)
		A
		B
		G
		P_0
		P_f
		Q
		R
		N
		chance_constraints_state = {}
		chance_constraints_control = {}
		mu_0 = []
		mu_f = []
		waypoints = {}
		objective_type = 'LQG'
		covariance_scaling = 1
		
		% Iteration parameters
		max_iters = 20  % Maximum number of iterations
		tol_opt = 1e-4  % Tolerance for objective function change between iterations
		tol_feas = 1e-4  % Tolerance for constraint feasibility violation
		verbose = true  % Print iteration information
		
		% Results
		prob  % Final FullCovarianceSteering problem instance
		iter_history = struct('iter', {}, 'P_ref', {}, 'Y_ref', {}, 'objective', {}, 'diagnostic', {}, 'change_objective', {}, 'max_violation', {})
	end
	
	methods
		function obj = FullCovarianceSteeringIterative(options)
			% Constructor
			% Accepts all parameters that FullCovarianceSteering accepts,
			% plus iteration parameters: max_iters, tol_change, verbose
			arguments
				options.A
				options.B
				options.G
				options.P_0
				options.P_f
				options.Q
				options.R
				options.N
				options.chance_constraints_state = {}
				options.chance_constraints_control = {}
				options.mu_0 = []
				options.mu_f = []
				options.waypoints = {}
				options.objective_type = 'LQG'
				options.covariance_scaling = 1
				options.max_iters = 20
				options.tol_opt = 1e-4
				options.tol_feas = 1e-4
				options.verbose = true
			end
			
			% Store problem parameters
			obj.A = options.A;
			obj.B = options.B;
			obj.G = options.G;
			obj.P_0 = options.P_0;
			obj.P_f = options.P_f;
			obj.Q = options.Q;
			obj.R = options.R;
			obj.N = options.N;
			obj.chance_constraints_state = options.chance_constraints_state;
			obj.chance_constraints_control = options.chance_constraints_control;
			obj.mu_0 = options.mu_0;
			obj.mu_f = options.mu_f;
			obj.waypoints = options.waypoints;
			obj.objective_type = options.objective_type;
			obj.covariance_scaling = options.covariance_scaling;
			
			% Store iteration parameters
			obj.max_iters = options.max_iters;
			obj.tol_opt = options.tol_opt;
			obj.tol_feas = options.tol_feas;
			obj.verbose = options.verbose;
		end
		
		function [diagnostic, prob] = solve(obj, sdp_settings)
			% Solve the problem iteratively
			% 
			% Inputs:
			%   sdp_settings - YALMIP solver settings (optional)
			% 
			% Outputs:
			%   diagnostic - Final diagnostic from YALMIP
			%   prob - Final FullCovarianceSteering problem instance
			
			if nargin < 2 || isempty(sdp_settings)
				sdp_settings = sdpsettings('verbose', 0, 'solver', 'mosek');
			end
			
			% Check if chance constraints are provided
			has_chance_constraints = ~isempty(obj.chance_constraints_state) || ~isempty(obj.chance_constraints_control);
			
			if ~has_chance_constraints
				% No chance constraints, solve once without iteration
				if obj.verbose
					fprintf('No chance constraints provided. Solving once without iteration.\n');
				end
				obj.prob = FullCovarianceSteering(...
					'A', obj.A, 'B', obj.B, 'G', obj.G, ...
					'P_0', obj.P_0, 'P_f', obj.P_f, ...
					'Q', obj.Q, 'R', obj.R, ...
					'N', obj.N, ...
					'mu_0', obj.mu_0, 'mu_f', obj.mu_f, ...
					'waypoints', obj.waypoints, ...
					'objective_type', obj.objective_type, ...
					'covariance_scaling', obj.covariance_scaling);
				
				diagnostic = obj.prob.solve(sdp_settings);
				
				if obj.verbose
					fprintf('  Exit status: %s\n', yalmiperror(diagnostic.problem));
				end
				
				if diagnostic.problem ~= 0 && diagnostic.problem ~= 4
					error('Solver failed: %s', yalmiperror(diagnostic.problem));
				end
				
				prob = obj.prob;
				return;
			end
			
			% Step 1: Solve without chance constraints to get initial reference
			if obj.verbose
				fprintf('=== Iteration 0: Solving without chance constraints ===\n');
			end
			
			prob_init = FullCovarianceSteering(...
				'A', obj.A, 'B', obj.B, 'G', obj.G, ...
				'P_0', obj.P_0, 'P_f', obj.P_f, ...
				'Q', obj.Q, 'R', obj.R, ...
				'N', obj.N, ...
				'mu_0', obj.mu_0, 'mu_f', obj.mu_f, ...
				'waypoints', obj.waypoints, ...
				'objective_type', obj.objective_type, ...
				'covariance_scaling', obj.covariance_scaling);
			
			diagnostic_init = prob_init.solve(sdp_settings);
			
			if obj.verbose
				fprintf('  Exit status: %s\n', yalmiperror(diagnostic_init.problem));
			end
			
			if diagnostic_init.problem ~= 0 && diagnostic_init.problem ~= 4
				error('Initial solve (without chance constraints) failed: %s', yalmiperror(diagnostic_init.problem));
			end
			
			% Extract reference from initial solution
			P_ref = prob_init.P;
			Y_ref = prob_init.P_u;
			if obj.verbose
				fprintf('  Initial objective: %.6f\n', prob_init.optimal_objective);
			end
			
			% Store initial iteration
			obj.iter_history(1).iter = 0;
			obj.iter_history(1).P_ref = P_ref;
			obj.iter_history(1).Y_ref = Y_ref;
			obj.iter_history(1).objective = prob_init.optimal_objective;
			obj.iter_history(1).diagnostic = diagnostic_init;
			obj.iter_history(1).change_objective = NaN;
			obj.iter_history(1).max_violation = NaN;
			
			% Step 2: Iteratively solve with chance constraints
			for iter = 1:obj.max_iters
				if obj.verbose
					fprintf('\n=== Iteration %d: Solving with chance constraints ===\n', iter);
				end
				
				% Create problem with current reference
				obj.prob = FullCovarianceSteering(...
					'A', obj.A, 'B', obj.B, 'G', obj.G, ...
					'P_0', obj.P_0, 'P_f', obj.P_f, ...
					'P_ref', P_ref, 'Y_ref', Y_ref, ...
					'Q', obj.Q, 'R', obj.R, ...
					'N', obj.N, ...
					'chance_constraints_state', obj.chance_constraints_state, ...
					'chance_constraints_control', obj.chance_constraints_control, ...
					'mu_0', obj.mu_0, 'mu_f', obj.mu_f, ...
					'waypoints', obj.waypoints, ...
					'objective_type', obj.objective_type, ...
					'covariance_scaling', obj.covariance_scaling);
				
				diagnostic = obj.prob.solve(sdp_settings);
				
				if obj.verbose
					fprintf('  Exit status: %s\n', yalmiperror(diagnostic.problem));
				end
				
				if diagnostic.problem ~= 0 && diagnostic.problem ~= 4
					error('Iteration %d failed: %s', iter, yalmiperror(diagnostic.problem));
				end
				
				% Compute objective function change
				obj_current = obj.prob.optimal_objective;
				obj_prev = obj.iter_history(iter).objective;
				if abs(obj_prev) > 1e-10
					change_objective = abs(obj_current - obj_prev) / abs(obj_prev);
				else
					change_objective = abs(obj_current - obj_prev);
				end
				
				% Compute constraint feasibility violations
				max_violation = 0;
				
				% Check state chance constraints
				if ~isempty(obj.chance_constraints_state) && ~isempty(obj.prob.mu)
					for i = 1:length(obj.chance_constraints_state)
						cc = obj.chance_constraints_state{i};
						if ~isfield(cc, 'type') || ~isfield(cc, 'p')
							continue;
						end
						
						if strcmp(cc.type, 'affine')
							alpha = cc.alpha(:);
							beta = cc.beta;
							p = cc.p;
							z = norminv(1 - p);
							
							% Determine which nodes to check
							if isfield(cc, 'nodes') && ~isempty(cc.nodes)
								nodes = cc.nodes;
							else
								nodes = 1:obj.N+1;
							end
							
							for k = nodes
								if k >= 1 && k <= obj.N+1
									mu_k = obj.prob.mu(:,k);
									P_k = obj.prob.P(:,:,k);
									% Ensure P_k is symmetric and compute alpha'*P_k*alpha
									P_k = (P_k + P_k') / 2;  % Symmetrize
									alpha_P_alpha = alpha' * P_k * alpha;
									
									% Check if P_k is positive semidefinite (alpha'*P_k*alpha >= 0)
									if alpha_P_alpha < 0
										% P_k is not PSD, set violation to a large value
										violation = inf;
									else
										% Actual constraint: alpha'*mu + z*sqrt(alpha'*P*alpha) <= beta
										constraint_val = alpha' * mu_k + z * sqrt(alpha_P_alpha);
										violation = max(0, constraint_val - beta);
									end
									max_violation = max(max_violation, violation);
								end
							end
						end
					end
				end
				
				% Check control chance constraints
				if ~isempty(obj.chance_constraints_control) && ~isempty(obj.prob.v)
					for i = 1:length(obj.chance_constraints_control)
						cc = obj.chance_constraints_control{i};
						if ~isfield(cc, 'type') || ~isfield(cc, 'p')
							continue;
						end
						
						if strcmp(cc.type, 'affine')
							alpha = cc.alpha(:);
							beta = cc.beta;
							p = cc.p;
							z = norminv(1 - p);
							
							% Determine which nodes to check
							if isfield(cc, 'nodes') && ~isempty(cc.nodes)
								nodes = cc.nodes;
							else
								nodes = 1:obj.N;
							end
							
							for k = nodes
								if k >= 1 && k <= obj.N
									v_k = obj.prob.v(:,k);
									Y_k = obj.prob.P_u(:,:,k);
									% Ensure Y_k is symmetric and compute alpha'*Y_k*alpha
									Y_k = (Y_k + Y_k') / 2;  % Symmetrize
									alpha_Y_alpha = alpha' * Y_k * alpha;
									
									% Check if Y_k is positive semidefinite (alpha'*Y_k*alpha >= 0)
									if alpha_Y_alpha < 0
										% Y_k is not PSD, set violation to a large value
										violation = inf;
									else
										% Actual constraint: alpha'*v + z*sqrt(alpha'*Y*alpha) <= beta
										constraint_val = alpha' * v_k + z * sqrt(alpha_Y_alpha);
										violation = max(0, constraint_val - beta);
									end
									max_violation = max(max_violation, violation);
								end
							end
	
						elseif strcmp(cc.type, 'norm')
							gamma = cc.gamma;
							p = cc.p;
							n = cc.n;
							q = sqrt(chi2inv(1 - p, n));
							
							% Determine which nodes to check
							if isfield(cc, 'nodes') && ~isempty(cc.nodes)
								nodes = cc.nodes;
							else
								nodes = 1:obj.N;
							end
							
							for k = nodes
								if k >= 1 && k <= obj.N
									v_k = obj.prob.v(:,k);
									Y_k = obj.prob.P_u(:,:,k);
									% Ensure Y_k is symmetric
									Y_k = (Y_k + Y_k') / 2;  % Symmetrize
									
									% Check if Y_k is positive semidefinite
									eig_Y = eig(Y_k);
									if any(eig_Y < -1e-10)  % Allow small numerical errors
										% Y_k is not PSD, set violation to a large value
										violation = inf;
									else
										% Actual constraint: norm(v) + q*sqrt(lambda_max(Y)) <= gamma
										lambda_max_Y = max(eig_Y);
										if lambda_max_Y < 0
											lambda_max_Y = 0;  % Clamp to 0 if negative due to numerical errors
										end
										constraint_val = norm(v_k, 2) + q * sqrt(lambda_max_Y);
										violation = max(0, constraint_val - gamma);
									end
									max_violation = max(max_violation, violation);
								end
							end
						end
					end
				end

				if obj.verbose
					fprintf('  Objective: %.6f\n', obj_current);
					fprintf('  Change in objective: %.6e\n', change_objective);
					fprintf('  Max constraint violation: %.6e\n', max_violation);
				end
				
				% Store iteration history
				obj.iter_history(iter+1).iter = iter;
				obj.iter_history(iter+1).P_ref = P_ref;
				obj.iter_history(iter+1).Y_ref = Y_ref;
				obj.iter_history(iter+1).objective = obj_current;
				obj.iter_history(iter+1).diagnostic = diagnostic;
				obj.iter_history(iter+1).change_objective = change_objective;
				obj.iter_history(iter+1).max_violation = max_violation;
				
				% Check convergence
				if change_objective < obj.tol_opt && max_violation < obj.tol_feas
					if obj.verbose
						fprintf('\nConverged after %d iterations.\n', iter);
					end
					break;
				end
				
				% Update reference for next iteration
				P_ref = obj.prob.P;
				Y_ref = obj.prob.P_u;
			end
			
			if iter >= obj.max_iters && obj.verbose
				fprintf('\nReached maximum iterations (%d).\n', obj.max_iters);
			end
			
			prob = obj.prob;
		end
	end
end