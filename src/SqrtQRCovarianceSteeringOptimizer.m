% Square root covariance steering optimizer with YALMIP optimizer capability
% This class inherits from SqrtQRCovarianceSteering and uses YALMIP's optimizer
% functionality to accelerate the algorithm by "pre-compiling" the optimization problem
classdef SqrtQRCovarianceSteeringOptimizer < SqrtQRCovarianceSteering

	properties 
		changing_params = struct();
		is_used_in_constraint_update = struct('S', true, 'L', true, 'mu', false, 'v', false);
	end

	methods
		function obj = SqrtQRCovarianceSteeringOptimizer(init_guess_struct, varargin)
			% Call parent constructor with all arguments
			% Note: setup_auxiliary_variables() will be called during initialize()
			% and will set up is_used_in_constraint_update for auxiliary variables
			obj@SqrtQRCovarianceSteering(init_guess_struct, varargin{:});
		end

		function val = get_use_optimizer(obj)
			val = true;
		end

		function is_aux = is_auxiliary_variable(obj, field_name)
			% Identify auxiliary variables used for SOC reformulation
			% These are variables introduced to reformulate quadratic objectives
			switch obj.objective_type
				case {'LQG', 'LQR', 'LQ'}
					auxiliary_fields = {'t_L', 't_S', 't_mu', 't_v'};
					is_aux = any(strcmp(field_name, auxiliary_fields));
				otherwise
					is_aux = false;
			end
		end

		function setup_auxiliary_variables(obj)
			% Override to also set up is_used_in_constraint_update for auxiliary variables
			% Call parent first to set impose_trust_region_struct
			setup_auxiliary_variables@SCPProblem(obj);
			
			% Set up is_used_in_constraint_update for auxiliary variables
			% Auxiliary variables are not used in constraint updates
			fields = fieldnames(obj.sdp_vars);
			for i = 1:length(fields)
				field = fields{i};
				if obj.is_auxiliary_variable(field)
					obj.is_used_in_constraint_update.(field) = false;
				end
			end
		end

		function initialize_auxiliary_variables(obj)
			% Initialize auxiliary variables based on primary variables
			% This sets initial values for t_L, t_S, t_mu, t_v based on L, S, mu, v
			switch obj.objective_type
				case {'LQG', 'LQR', 'LQ'}
					% Initialize auxiliary variables if they exist in init_guess_struct
					if isfield(obj.init_guess_struct, 'L') && isfield(obj.init_guess_struct, 'S')
						% Initialize t_L and t_S based on initial guess
						if ~isfield(obj.init_guess_struct, 't_L')
							R_sqrt_upper = chol(obj.R);
							obj.init_guess_struct.t_L = zeros(1, obj.N);
							for k = 1:obj.N
								if size(obj.init_guess_struct.L, 3) >= k
									L_k = obj.init_guess_struct.L(:,:,k);
									obj.init_guess_struct.t_L(k) = norm(R_sqrt_upper * L_k, 'fro')^2;
								end
							end
						end
						
						if ~isfield(obj.init_guess_struct, 't_S')
							Q_sqrt_upper = chol(obj.Q);
							obj.init_guess_struct.t_S = zeros(1, obj.N);
							for k = 1:obj.N
								if size(obj.init_guess_struct.S, 3) >= k
									S_k = obj.init_guess_struct.S(:,:,k);
									obj.init_guess_struct.t_S(k) = norm(Q_sqrt_upper * S_k, 'fro')^2;
								end
							end
						end
					end
					
					% Initialize t_mu and t_v if mean variables exist
					if isfield(obj.init_guess_struct, 'mu') && isfield(obj.init_guess_struct, 'v') ...
							&& ~isempty(obj.init_guess_struct.mu) && ~isempty(obj.init_guess_struct.v)
						if ~isfield(obj.init_guess_struct, 't_mu')
							Q_sqrt_upper = chol(obj.Q);
							obj.init_guess_struct.t_mu = zeros(1, obj.N);
							for k = 1:obj.N
								if size(obj.init_guess_struct.mu, 2) >= k
									mu_k = obj.init_guess_struct.mu(:,k);
									obj.init_guess_struct.t_mu(k) = norm(Q_sqrt_upper * mu_k)^2;
								end
							end
						end
						
						if ~isfield(obj.init_guess_struct, 't_v')
							R_sqrt_upper = chol(obj.R);
							obj.init_guess_struct.t_v = zeros(1, obj.N);
							for k = 1:obj.N
								if size(obj.init_guess_struct.v, 2) >= k
									v_k = obj.init_guess_struct.v(:,k);
									obj.init_guess_struct.t_v(k) = norm(R_sqrt_upper * v_k)^2;
								end
							end
						end
					end
			end
		end

		function vars = define_vars(obj)
			% Call parent to get base variables
			vars = define_vars@SqrtQRCovarianceSteering(obj);
			
			switch obj.objective_type
				case {'LQG', 'LQR', 'LQ'}
					
				% Add scalar variables for squared norms (to use SOC formulation)
				% These are needed to reformulate quadratic objectives for YALMIP optimizer pattern
				vars.t_L = sdpvar(1, obj.N, 'full');
				vars.t_S = sdpvar(1, obj.N, 'full');
				
				% Scalar variables for mean state and control cost (if mean variables exist)
				if isfield(vars, 'mu') && isfield(vars, 'v') && ~isempty(vars.mu) && ~isempty(vars.v)
					vars.t_mu = sdpvar(1, obj.N, 'full');
					vars.t_v = sdpvar(1, obj.N, 'full');
				end
			end
		end

		function J0 = objective(obj, vars)
			% Reformulate objective using auxiliary variables with SOC constraints
			% to avoid quadratic terms that become quadratic after fixing parameters
			J0 = 0;
			
			switch obj.objective_type
				case {'LQG', 'LQR', 'LQ'}
					% Use scalar variables for squared norms (SOC constraints added in convex_ineq)
					J0 = J0 + sum(vars.t_L) + sum(vars.t_S);
					
					% Mean state and control cost (if mean variables are defined)
					if isfield(vars, 'mu') && isfield(vars, 'v') && ~isempty(vars.mu) && ~isempty(vars.v)
						J0 = J0 + sum(vars.t_mu) + sum(vars.t_v);
					end
				case 'DV99'
					% For DV99, use the parent implementation (uses norm, which is fine)
					J0 = objective@SqrtQRCovarianceSteering(obj, vars);
				otherwise
					error('Unknown objective type');
			end
		end

		function constraints = convex_ineq(obj, vars)
			% Call parent to get base constraints
			constraints = convex_ineq@SqrtQRCovarianceSteering(obj, vars);
			
			% Add second-order cone constraints for squared norms
			% This reformulates ||A*x||^2 <= t as a SOC constraint: norm([2*A*x; t-1]) <= t+1
			% See Problem 4.26 in Boyd's book
			
			switch obj.objective_type
				case {'LQG', 'LQR', 'LQ'}
				R_sqrt_upper = chol(obj.R);
				Q_sqrt_upper = chol(obj.Q);
				
				% cone_matrix = [];
                % 
				% for k = 1:obj.N
				% 	cone_matrix(:,k) = 2 * reshape(R_sqrt_upper * vars.L(:,:,k), [], 1);
				% end
				% % cone_matrix(9,:) = vars.t_L - 1;
				% cone_matrix = [vars.t_L + 1; cone_matrix];
				% constraints = [constraints; cone(cone_matrix)];

				% SOC constraints for control cost: ||R_sqrt_upper * L||^2 <= t_L
				for k = 1:obj.N
					constraints = [constraints;
						% norm([2 * reshape(R_sqrt_upper * vars.L(:,:,k), [], 1); vars.t_L(k) - 1]) <= vars.t_L(k) + 1
						cone([2 * reshape(R_sqrt_upper * vars.L(:,:,k), [], 1); vars.t_L(k) - 1], vars.t_L(k) + 1)
					];
				end
				
				% SOC constraints for state cost: ||Q_sqrt_upper * S||^2 <= t_S
				for k = 1:obj.N
					constraints = [constraints;
						% norm([2 * reshape(Q_sqrt_upper * vars.S(:,:,k), [], 1); vars.t_S(k) - 1]) <= vars.t_S(k) + 1
						cone([2 * reshape(Q_sqrt_upper * vars.S(:,:,k), [], 1); vars.t_S(k) - 1], vars.t_S(k) + 1)
					];
				end

				constraints = [constraints;
					vars.t_L >= 0
					vars.t_S >= 0
				];
				
				% SOC constraints for mean state and control cost
				if isfield(vars, 'mu') && isfield(vars, 'v') && ~isempty(vars.mu) && ~isempty(vars.v)
					constraints = [constraints;
						cone([vars.t_mu + 1; 2 * Q_sqrt_upper * vars.mu(:,1:obj.N); vars.t_mu - 1])
						cone([vars.t_v + 1; 2 * R_sqrt_upper * vars.v; vars.t_v - 1])
					];
					% for k = 1:obj.N
					% 	% Mean state: ||Q_sqrt_upper * mu||^2 <= t_mu
					% 	constraints = [constraints;
					% 		% norm([2 * Q_sqrt_upper * vars.mu(:,k); vars.t_mu(k) - 1]) <= vars.t_mu(k) + 1
					% 		cone([2 * Q_sqrt_upper * vars.mu(:,k); vars.t_mu(k) - 1], vars.t_mu(k) + 1)
					% 	];
						
					% 	% Mean control: ||R_sqrt * v||^2 <= t_v
					% 	constraints = [constraints;
					% 		% norm([2 * R_sqrt_upper * vars.v(:,k); vars.t_v(k) - 1]) <= vars.t_v(k) + 1
					% 		cone([2 * R_sqrt_upper * vars.v(:,k); vars.t_v(k) - 1], vars.t_v(k) + 1)
					% 	];
					% end

					constraints = [constraints;
						vars.t_mu >= 0
						vars.t_v >= 0
					];
				end
			end
		end

		function set_changing_parameters_sdp(obj)
			% Define SDP variables for parameters that change between iterations
			% These correspond to the linearization points for the QR decomposition
			
			% Parameters for QR decomposition linearization
			obj.changing_params.Q_ref = sdpvar(obj.nx+obj.nw, obj.nx, obj.N, 'full');
			obj.changing_params.R_ref = sdpvar(obj.nx, obj.nx, obj.N, 'full');
			obj.changing_params.R_ref_inv = sdpvar(obj.nx, obj.nx, obj.N, 'full');
			for k = 1:obj.N
				obj.changing_params.R_ref(:,:,k) = triu(obj.changing_params.R_ref(:,:,k));
				obj.changing_params.R_ref_inv(:,:,k) = triu(obj.changing_params.R_ref_inv(:,:,k));
			end
		end

		function p = get_changing_parameters(obj, ref_vars)
			% Numerical update of the parameters that appear in the constraints
			% This function computes the reference values for linearization
			
			p = struct('Q_ref', zeros(obj.nx+obj.nw, obj.nx, obj.N), 'R_ref', zeros(obj.nx, obj.nx, obj.N), 'R_ref_inv', zeros(obj.nx, obj.nx, obj.N));
			
			for k = 1:obj.N
				A_k = obj.A_sys(:,:,k);
				B_k = obj.B_sys(:,:,k);
				G_k = obj.G_sys(:,:,k);
				S_ref_k = ref_vars.S(:,:,k);
				L_ref_k = ref_vars.L(:,:,k);

				X_ref_k = [A_k * S_ref_k + B_k * L_ref_k, G_k];
				[p.Q_ref(:,:,k), p.R_ref(:,:,k)] = obj.economy_qr_with_positive_diagonal(X_ref_k');
				p.R_ref_inv(:,:,k) = inv(p.R_ref(:,:,k));
			end
		end

		function constraints = get_changing_constraints(obj, vars, ref_vars)
			% Get the changing constraints using YALMIP optimizer
			obj.set_changing_parameters_sdp();
			p = obj.changing_params;

			% Linearized QR decomposition constraints
			% The constraints are built in the same order as noncvx_eq_relaxed
			constraintLHS = [];
			
			for k = 1:obj.N
				% Get system matrices for time step k
				A_k = obj.A_sys(:,:,k);
				B_k = obj.B_sys(:,:,k);
				G_k = obj.G_sys(:,:,k);
				
				% Current and next covariance square roots
				S_k = vars.S(:,:,k);
				S_kp1 = vars.S(:,:,k+1);
				L_k = vars.L(:,:,k);
				
				% Reference values for linearization
				S_k_ref = ref_vars.S(:,:,k);
				L_k_ref = ref_vars.L(:,:,k);
				
				% Build the matrix for QR decomposition at reference point
				X_k_ref = [A_k * S_k_ref + B_k * L_k_ref, G_k];
				dX = [A_k * (S_k - S_k_ref) + B_k * (L_k - L_k_ref), zeros(size(G_k))];
				
				% Compute QR derivative using the d_QR function
				[dR, ~, R_X_k_ref] = obj.d_QR(X_k_ref', dX', p.Q_ref(:,:,k), p.R_ref(:,:,k), p.R_ref_inv(:,:,k));

				% Add linearized constraint (same structure as noncvx_eq_relaxed)
				constraintLHS = [constraintLHS
					obj.vec_tril(S_kp1 - R_X_k_ref' - dR')
				];
			end
			
			% The constraint is that the linearized constraint equals the slack variable
			constraints = [constraintLHS == obj.slack_noncvx_eq];
		end

	end

end
