classdef FullCovarianceSteering < CovarianceSteeringBase
	% Full Covariance Steering class
	% 
	% This class implements covariance steering using full covariance matrices
	% as optimization variables.
	%
	% Usage:

	%   cs = FullCovarianceSteering('A', A, 'B', B, 'G', G, 'P_0', P_0, 'P_f', P_f, ...
	%                               'Q', Q, 'R', R, 'nx', nx, 'nu', nu, 'N', N);
	%   diagnostic = cs.solve_problem();
	
	methods
		function obj = FullCovarianceSteering(varargin)
			obj = obj@CovarianceSteeringBase(varargin{:});
		end

		function set_sdpvars(obj)
			% Define YALMIP sdpvar variables
			nx = obj.nx;
			nu = obj.nu;
			N = obj.N;

			obj.sdpvars.P = sdpvar(nx, nx, N+1); % State covariance
			obj.sdpvars.U = sdpvar(nu, nx, N); % Control gain
			obj.sdpvars.Y = sdpvar(nu, nu, N); % Auxiliary variable for input covariance
			
			% Mean variables (if mean constraints are specified)
			if ~isempty(obj.chance_constraints_state) || ~isempty(obj.chance_constraints_control) ...
					|| ~isempty(obj.mu_0) || ~isempty(obj.mu_f) || ~isempty(obj.waypoints)
				obj.sdpvars.mu = sdpvar(nx, N+1, 'full'); % State mean trajectory
				obj.sdpvars.v = sdpvar(nu, N, 'full'); % Control mean trajectory
			else
				obj.sdpvars.mu = [];
				obj.sdpvars.v = [];
			end
			
			obj.sdpvars.J = [];
			obj.sdpvars.constraints = [];
		end

		function set_objective(obj)
			switch obj.objective_type
				case {'LQG', 'LQR', 'LQ'}
					% Covariance part of objective (always included)
					J_cov = 0;
					for k = 1:obj.N
						J_cov = J_cov + trace(obj.Q * obj.sdpvars.P(:,:,k)) ...
							+ trace(obj.R * obj.sdpvars.Y(:,:,k));
					end
					
					% Mean part of objective (if mean variables are defined)
					if ~isempty(obj.sdpvars.mu) && ~isempty(obj.sdpvars.v)
						J_mean = 0;
						for k = 1:obj.N
							J_mean = J_mean + obj.sdpvars.mu(:,k)' * obj.Q * obj.sdpvars.mu(:,k) ...
								+ obj.sdpvars.v(:,k)' * obj.R * obj.sdpvars.v(:,k);
						end
						% Add terminal state mean cost
						J_mean = J_mean + obj.sdpvars.mu(:,obj.N+1)' * obj.Q * obj.sdpvars.mu(:,obj.N+1);
						obj.sdpvars.J = J_mean + J_cov;
					else
						obj.sdpvars.J = J_cov;
					end
					
				case 'DV99'
					% DV99 objective: minimize 99th percentile of control norm
					% Formulation: sum_k [norm(v_k) + q * sqrt(lambda_max(Y_k))]
					% where q = sqrt(chi2inv(0.99, nu))
					% Since sqrt(lambda_max(Y_k)) is not convex, we linearize around Y_ref
					% Similar to norm control chance constraints

					% add a small quadratic term to the objective to ensure losslessness
					epsilon_lossless = 1E-3;
					if isempty(obj.Y_ref)
						error('FullCovarianceSteering: Y_ref must be provided for DV99 objective.');
					end
					
					if isempty(obj.sdpvars.v)
						error('FullCovarianceSteering: Mean control variables (v) must be defined for DV99 objective.');
					end
					
					q = sqrt(chi2inv(0.99, obj.nu));
					obj.sdpvars.J = 0;
					
					for k = 1:obj.N
						obj.sdpvars.J = obj.sdpvars.J + epsilon_lossless * trace(obj.sdpvars.Y(:,:,k)) / obj.covariance_scaling;

						% Get reference for linearization
						if size(obj.Y_ref, 3) == 1
							Y_ref_k = obj.Y_ref;
						else
							Y_ref_k = obj.Y_ref(:,:,k);
						end
						Y_ref_k = Y_ref_k * obj.covariance_scaling;
						sqrt_lambda_max_ref = sqrt(lambda_max(Y_ref_k));
						
						% Linearized objective: norm(v_k) + q * [sqrt_lambda_max_ref + lambda_max(Y_k) / (2 * sqrt_lambda_max_ref)]
						obj.sdpvars.J = obj.sdpvars.J + norm(obj.sdpvars.v(:,k), 2) ...
							+ sqrt(obj.covariance_scaling) * q * sqrt_lambda_max_ref ...
							+ sqrt(obj.covariance_scaling) * q * lambda_max(obj.sdpvars.Y(:,:,k)) / (2 * sqrt_lambda_max_ref);
					end
					
				otherwise
					error('FullCovarianceSteering: Unknown objective type ''%s''.', obj.objective_type);
			end
		end

		function set_constraints(obj)
			% Set all constraints at once
			constraints = [
				obj.get_boundary_constraints()
				obj.get_dynamics_constraints()
				obj.get_mean_dynamics_constraints()
				obj.get_LMI_constraints(1:obj.N)
				obj.get_positive_semidefinite_constraints()
				obj.get_mean_boundary_constraints()
				obj.get_chance_constraints_state()
				obj.get_chance_constraints_control()
			];

			obj.sdpvars.constraints = constraints;
		end

		function constraints = get_boundary_constraints(obj)
			constraints = [];
			constraints = [constraints
				[obj.sdpvars.P(:,:,1) == obj.P_0 .* obj.covariance_scaling]:'initial_covariance'
				[obj.P_f .* obj.covariance_scaling >= obj.sdpvars.P(:,:,obj.N+1)]:'final_covariance'
			];
		end
		
		function constraints = get_mean_dynamics_constraints(obj)
			% Mean dynamics: mu_{k+1} = A_k * mu_k + B_k * v_k
			constraints = [];
			if ~isempty(obj.sdpvars.mu) && ~isempty(obj.sdpvars.v)
				for k = 1:obj.N
					constraints = [constraints
						obj.sdpvars.mu(:,k+1) == obj.A(:,:,k) * obj.sdpvars.mu(:,k) + obj.B(:,:,k) * obj.sdpvars.v(:,k)
					];
				end
			end
		end
		
		function constraints = get_mean_boundary_constraints(obj)
			% Mean boundary conditions and waypoints
			constraints = [];
			if ~isempty(obj.sdpvars.mu)
				if ~isempty(obj.mu_0)
					constraints = [constraints
						obj.sdpvars.mu(:,1) == obj.mu_0(:)
					];
				end
				if ~isempty(obj.mu_f)
					constraints = [constraints
						obj.sdpvars.mu(:,obj.N+1) == obj.mu_f(:)
					];
				end
				% Waypoint constraints (intermediate nodes)
				% NaN values in wp.mu indicate unconstrained components
				if ~isempty(obj.waypoints)
					for i = 1:length(obj.waypoints)
						wp = obj.waypoints{i};
						if isfield(wp, 'node') && isfield(wp, 'mu') && wp.node >= 1 && wp.node <= obj.N+1
							mu_wp = wp.mu(:);
							% Find non-NaN components to constrain
							idx_constrained = ~isnan(mu_wp);
							if any(idx_constrained)
								constraints = [constraints
									obj.sdpvars.mu(idx_constrained, wp.node) == mu_wp(idx_constrained)
								];
							end
						end
					end
				end
			end
		end
		
		function constraints = get_chance_constraints_state(obj)
			% State chance constraints for full covariance method
			constraints = [];
			if ~isempty(obj.chance_constraints_state) && ~isempty(obj.sdpvars.mu)
				for i = 1:length(obj.chance_constraints_state)
					cc = obj.chance_constraints_state{i};
					if ~isfield(cc, 'type') || ~isfield(cc, 'p')
						continue;
					end
					
					% Determine which nodes to apply constraint to
					if isfield(cc, 'nodes') && ~isempty(cc.nodes)
						nodes = cc.nodes;
					else
						nodes = 1:obj.N+1; % Apply to all nodes if not specified
					end
					
					switch cc.type
						case 'affine'
							% Affine chance constraint: P(alpha'*x <= beta) >= 1-p
							% Using linearized formulation around reference covariance P_ref:
							% z * (1/(2*sqrt(alpha'*P_ref*alpha))) * alpha'*P_k*alpha + alpha'*mu_k 
							%   - (beta - z*(1/2)*sqrt(alpha'*P_ref*alpha)) <= 0
							% This is the first-order Taylor expansion of sqrt(alpha'*P_k*alpha) around P_ref
							alpha = cc.alpha(:); % ensure column vector
							beta = cc.beta;
							p = cc.p;
							z = norminv(1 - p);

							if isempty(obj.P_ref)
								error('FullCovarianceSteering: P_ref must be provided for affine chance constraints.');
							end
							
							for k = nodes
								if k >= 1 && k <= obj.N+1
		
									if size(obj.P_ref, 3) == 1
										P_ref_k = obj.P_ref;
									else
										P_ref_k = obj.P_ref(:,:,k);
									end
									sqrt_ref = sqrt(alpha' * P_ref_k * alpha);
									constraints = [constraints
										z / (2 * sqrt_ref) * (alpha' * obj.sdpvars.P(:,:,k) * alpha) ...
											+ alpha' * obj.sdpvars.mu(:,k) - beta + z * sqrt_ref / 2 <= 0
									];
								end
							end
					
						% case 'norm'
						% 	% Norm chance constraint: P(||x||_2 <= gamma) >= 1-p
						% 	gamma = cc.gamma;
						% 	p = cc.p;
						% 	n = cc.n;
						% 	q = sqrt(chi2inv(1 - p, n));
							
						% 	error('FullCovarianceSteering: Norm chance constraints are not supported for full covariance method.');
						otherwise
							error('FullCovarianceSteering: Unsupported chance constraint type.');
					end
				end
			end
		end
		
		function constraints = get_chance_constraints_control(obj)
			% Control chance constraints for full covariance method
			constraints = [];
			if ~isempty(obj.chance_constraints_control) && ~isempty(obj.sdpvars.v)
				for i = 1:length(obj.chance_constraints_control)
					cc = obj.chance_constraints_control{i};
					if ~isfield(cc, 'type') || ~isfield(cc, 'p')
						continue;
					end
					
					% Determine which nodes to apply constraint to
					if isfield(cc, 'nodes') && ~isempty(cc.nodes)
						nodes = cc.nodes;
					else
						nodes = 1:obj.N; % Apply to all control time steps if not specified
					end
									
					if isempty(obj.Y_ref)
						error('FullCovarianceSteering: Y_ref must be provided for affine control chance constraints.');
					end

					switch cc.type
						case 'affine'
							% Affine chance constraint: P(alpha'*u <= beta) >= 1-p
							% Using linearized formulation around reference control covariance Y_ref:
							% z * (1/(2*sqrt(alpha'*Y_ref*alpha))) * alpha'*Y_k*alpha + alpha'*v_k 
							%   - (beta - z*(1/2)*sqrt(alpha'*Y_ref*alpha)) <= 0
							% This is the first-order Taylor expansion of sqrt(alpha'*Y_k*alpha) around Y_ref
							alpha = cc.alpha(:); % ensure column vector
							beta = cc.beta;
							p = cc.p;
							z = norminv(1 - p);
							
							for k = nodes
								if k >= 1 && k <= obj.N
									if size(obj.Y_ref, 3) == 1
										Y_ref_k = obj.Y_ref;
									else
										Y_ref_k = obj.Y_ref(:,:,k);
									end
									sqrt_ref = sqrt(alpha' * Y_ref_k * alpha);
									constraints = [constraints
										z / (2 * sqrt_ref) * (alpha' * obj.sdpvars.Y(:,:,k) * alpha) ...
											+ alpha' * obj.sdpvars.v(:,k) - beta + z * sqrt_ref / 2 <= 0
									];
								end
							end
							
						case 'norm'
							% Norm chance constraint: P(||u||_2 <= gamma) >= 1-p
							% Formulation: norm(v_k, 2) + q * sqrt(lambda_max(Y_k)) <= gamma
							% where q = sqrt(chi2inv(1-p, n))
							gamma = cc.gamma;
							p = cc.p;
							n = cc.n;
							q = sqrt(chi2inv(1 - p, n));
							
							for k = nodes
								if k >= 1 && k <= obj.N
									% Conservative approximation: norm(v_k) + q * sqrt(lambda_max(Y_k)) <= gamma
									if size(obj.Y_ref, 3) == 1
										Y_ref_k = obj.Y_ref;
									else
										Y_ref_k = obj.Y_ref(:,:,k);
									end
									sqrt_lambda_max_ref = sqrt(lambda_max(Y_ref_k));
									constraints = [constraints
										norm(obj.sdpvars.v(:,k), 2) + q * sqrt_lambda_max_ref + q * lambda_max(obj.sdpvars.Y(:,:,k)) / (2 * sqrt_lambda_max_ref) - gamma <= 0
									];
								end
							end
						otherwise
							error('FullCovarianceSteering: Unsupported control chance constraint type.');
					end
				end
			end
		end

		function constraints = get_dynamics_constraints(obj)
			A = obj.A;
			B = obj.B;
			G = obj.G;

			constraints = [];
			for k = 1:obj.N
				constraints = [constraints
					obj.sdpvars.P(:,:,k+1) == ...
						(A(:,:,k) * obj.sdpvars.P(:,:,k) * A(:,:,k)' ...
						+ A(:,:,k) * obj.sdpvars.U(:,:,k)' * B(:,:,k)' ...
						+ B(:,:,k) * obj.sdpvars.U(:,:,k) * A(:,:,k)' ...
						+ B(:,:,k) * obj.sdpvars.Y(:,:,k) * B(:,:,k)' ...
						+ (G(:,:,k) * G(:,:,k)') .* obj.covariance_scaling)
				];
			end
			constraints = constraints:'dynamics';
		end

		function constraints = get_LMI_constraints(obj, node_indices)
			constraints = [];
			for k = node_indices
				constraints = [constraints
					[ obj.sdpvars.P(:,:,k) , obj.sdpvars.U(:,:,k)' ;
						obj.sdpvars.U(:,:,k) , obj.sdpvars.Y(:,:,k) ] >= 0
				];
			end
			constraints = constraints:'LMI';
		end

		function constraints = get_positive_semidefinite_constraints(obj)
			constraints = [];
			for k = 1:obj.N
				constraints = [constraints
					[obj.sdpvars.Y(:,:,k) >= 0]:'Y_positive_semidefinite'
				];
			end
		end

		function diagnostic = solve_problem(obj, options)
			% Override solve_problem to add losslessness check for full covariance
			arguments
				obj
				options.verbose = 0
				options.solver = 'mosek'
				options.savesolveroutput = true
			end
			
			diagnostic = solve_problem@CovarianceSteeringBase(obj, options);
			
			% Additional losslessness check for full covariance method
			if diagnostic.problem == 0 || diagnostic.problem == 4
				if options.verbose
					fprintf('Optimization with full_covariance successful. Objective value: %f\n', obj.optimal_objective);
				end
				if diagnostic.problem == 0
					obj.check_lossless(obj.P, obj.K, obj.A, obj.B, obj.G, obj.N);
				elseif diagnostic.problem == 4
					fprintf("Numerical problems. Objective value: %f.\n", obj.optimal_objective);
					obj.check_lossless(obj.P, obj.K, obj.A, obj.B, obj.G, obj.N);
				end
			end
		end

		function set_feedback_gains(obj)
			% Compute feedback gains from optimal U and P
			nx = obj.nx;
			nu = obj.nu;
			N = obj.N;
			K = zeros(nu, nx, N);
			
			% Store state and control covariances
			obj.P = value(obj.sdpvars.P) / obj.covariance_scaling;
			obj.P_u = value(obj.sdpvars.Y) / obj.covariance_scaling;
			
			for k = 1:N
				P_k = obj.P(:,:,k);
				U_k = value(obj.sdpvars.U(:,:,k)) / obj.covariance_scaling;
				K(:,:,k) = U_k / P_k; % Feedback gain
			end
			obj.K = K;
			
			% Store mean trajectories if they exist
			if ~isempty(obj.sdpvars.mu)
				obj.mu = value(obj.sdpvars.mu);
			else
				obj.mu = [];
			end
			
			if ~isempty(obj.sdpvars.v)
				obj.v = value(obj.sdpvars.v);
			else
				obj.v = [];
			end
		end

		function [is_lossless, worst_loss] = check_lossless(obj, tol, verbose)
			% Check losslessness of the solution for full covariance method
			arguments
				obj
				tol = 1E-4
				verbose = false
			end

			if isempty(obj.P)
				error('FullCovarianceSteering: P not set. Solve the problem first.');
			end
			K = obj.K;
			A = obj.A;
			B = obj.B;
			G = obj.G;
			N = obj.N;

			worst_loss = 0;
			is_lossless = true;
			for k = 1:N
				P_next = (A(:,:,k) + B(:,:,k)*K(:,:,k)) * obj.P(:,:,k) * (A(:,:,k) + B(:,:,k)*K(:,:,k))' + G(:,:,k)*G(:,:,k)';
				% Add control-dependent noise term if available
	
				loss = norm(P_next - obj.P(:,:,k+1), 'fro') / norm(obj.P(:,:,k+1), 'fro');
				if loss > tol
					is_lossless = false;
					worst_loss = max(worst_loss, loss);
				end
			end

			if is_lossless && verbose
				fprintf('Losslessness verified within tolerance %g.\n', tol);
			else
				fprintf('Losslessness NOT verified. Worst relative loss: %g\n', worst_loss);
			end
		end

		function [is_SDP] = control_covariance_is_SDP(obj, verbose)
			arguments
				obj
				verbose = false
			end
			is_SDP = true;
			for k = 1:obj.N
				[~, flag] = chol(obj.P_u(:,:,k));
				is_SDP = is_SDP && flag == 0;
			end
			if verbose
				if is_SDP
					fprintf('Control covariance is SDP.\n');
				else
					fprintf('Control covariance is not SDP.\n');
				end
			end
		end
	end
end