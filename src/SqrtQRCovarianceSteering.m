% Square root covariance steering with QR decomposition-based covariance propagation
% This class inherits from SCPProblem 
classdef SqrtQRCovarianceSteering < SCPProblem

	properties
		init_guess_struct
		impose_trust_region_struct = struct('S', true, 'L', true, 'mu', false, 'v', false);
	end

	properties
		N % Horizon length
		nx % State dimension
		nu % Control dimension
		nw % Process noise dimension
		A_sys % System dynamics matrices
		B_sys % Control input matrices
		G_sys % Process noise matrices
		P_0 % Initial state covariance
		P_f % Terminal state covariance
		Q % State cost matrix
		R % Control cost matrix
		P % state covariance, set after solving
		P_u % control covariance, set after solving
		K % feedback gain, set after solving
		objective_type = 'LQG'; % 'LQG' or 'DV99'
		% Chance constraint support
		% State constraints: cell array of structs with fields:
		%   - type: 'affine' or 'norm'
		%   - For affine: alpha (vector), beta (scalar), p (violation prob), nodes (optional, default: all)
		%   - For norm: gamma (scalar), p (violation prob), n (dimension), nodes (optional, default: all)
		chance_constraints_state = {}
		chance_constraints_control = {} % Same format as state constraints
		circular_obstacles = {} % cell array of circular obstacles with fields: 'center' (nx x 1 vector), 'radius' (scalar), 'p' (violation prob)
		mu_0 = []  % initial mean (nx x 1) or empty -> assumed zero
		mu_f = []  % terminal mean (nx x 1) or empty -> assumed zero
		waypoints = {}  % cell array of waypoint structs with fields: 'node' (scalar, 1:N+1) and 'mu' (nx x 1 vector)
		mu % state mean trajectory, set after solving
		v  % control mean trajectory, set after solving

	end

	methods
		function obj = SqrtQRCovarianceSteering(init_guess_struct, options)
			arguments
				init_guess_struct struct
				options.N
				options.A_sys
				options.B_sys
				options.G_sys
				options.P_0
				options.P_f
				options.Q = []
				options.R = []
				options.objective_type = 'LQG';
				options.chance_constraints_state = {}
				options.chance_constraints_control = {}
				options.mu_0 = []
				options.mu_f = []
				options.waypoints = {}
				options.D = [] % Trust region scaling matrix. If empty, defaults to scalar 1
				options.circular_obstacles = {}
			end
			obj@SCPProblem();
			obj.init_guess_struct = init_guess_struct;

			obj.N = options.N;
			obj.nx = size(options.A_sys, 1);
			obj.nu = size(options.B_sys, 2);
			obj.nw = size(options.G_sys, 2);
			obj.A_sys = options.A_sys;
			obj.B_sys = options.B_sys;
			obj.G_sys = options.G_sys;
			obj.P_0 = options.P_0;
			obj.P_f = options.P_f;
			obj.Q = options.Q;
			obj.R = options.R;
			obj.objective_type = options.objective_type;
			% optional chance constraints and endpoint means
			obj.chance_constraints_state = options.chance_constraints_state;
			obj.chance_constraints_control = options.chance_constraints_control;
			obj.mu_0 = options.mu_0;
			obj.mu_f = options.mu_f;
			obj.waypoints = options.waypoints;
			% Trust region scaling matrix
			obj.D = options.D;
			obj.circular_obstacles = options.circular_obstacles;

			yalmip('clear');
			obj.initialize();

		end

		function vars = define_vars(obj)
			vars.S = sdpvar(obj.nx, obj.nx, obj.N+1);
			for k = 1:obj.N+1
				vars.S(:,:,k) = tril(vars.S(:,:,k));
			end
			vars.L = sdpvar(obj.nu, obj.nx, obj.N);
			
			% Mean variables (only if needed for chance constraints, boundary conditions, or waypoints)
			if ~isempty(obj.chance_constraints_state) || ~isempty(obj.chance_constraints_control) ...
					|| ~isempty(obj.mu_0) || ~isempty(obj.mu_f) || ~isempty(obj.waypoints)
				vars.mu = sdpvar(obj.nx, obj.N+1);
				vars.v = sdpvar(obj.nu, obj.N);
			else
				vars.mu = [];
				vars.v = [];
			end
		end

		function J0 = objective(obj, vars)
			% Minimize control effort and state covariance
			J0 = 0;
			
			switch obj.objective_type
				case {'LQG', 'LQR', 'LQ'}

					% Add control cost
					for k = 1:obj.N
						J0 = J0 + trace(vars.L(:,:,k) * vars.L(:,:,k)' * obj.R);
					end
					
					% Add state covariance cost
					for k = 1:obj.N
						J0 = J0 + trace(vars.S(:,:,k) * vars.S(:,:,k)' * obj.Q);
					end

					% Add mean state and control cost (if mean variables are defined)
					if isfield(vars, 'mu') && isfield(vars, 'v') && ~isempty(vars.mu) && ~isempty(vars.v)
						for k = 1:obj.N
							J0 = J0 + vars.mu(:,k)' * obj.Q * vars.mu(:,k) + vars.v(:,k)' * obj.R * vars.v(:,k);
						end
					end
				case 'DV99'
					% Add control cost only
					q = sqrt(chi2inv(0.99, obj.nu));
					for k = 1:obj.N
						if isfield(vars, 'v') && ~isempty(vars.v)
							J0 = J0 + norm(vars.v(:,k)) + q * norm(vars.L(:,:,k), 2);
						else
							J0 = J0 + q * norm(vars.L(:,:,k), 2);
						end
					end
				otherwise
					error('Unknown objective type');
			end
		end

		function constraints = convex_eq(obj, vars)
			constraints = [
				[vars.S(:,:,1) == chol(obj.P_0, 'lower')]:'Initial Covariance'
			];

			% Mean dynamics constraints (only if mean variables are defined)
			if isfield(vars, 'mu') && isfield(vars, 'v') && ~isempty(vars.mu) && ~isempty(vars.v)
				% initial mean if provided, otherwise zero
				if ~isempty(obj.mu_0)
					constraints = [constraints; [vars.mu(:,1) == obj.mu_0]:'Initial Mean'];
				end
				
				% Mean dynamics: mu_{k+1} = A_k * mu_k + B_k * v_k
				for k = 1:obj.N
					A_k = obj.A_sys(:,:,k);
					B_k = obj.B_sys(:,:,k);
					constraints = [constraints; [vars.mu(:,k+1) == A_k * vars.mu(:,k) + B_k * vars.v(:,k)]:'Mean Dynamics'];
				end
				
				% enforce terminal mean if provided in mu_f
				if ~isempty(obj.mu_f)
					constraints = [constraints; [vars.mu(:,obj.N+1) == obj.mu_f]:'Terminal Mean'];
				end
				
				% enforce waypoint mean constraints (intermediate nodes)
				% NaN values in wp.mu indicate unconstrained components
				if ~isempty(obj.waypoints)
					for i = 1:length(obj.waypoints)
						wp = obj.waypoints{i};
						if isfield(wp, 'node') && isfield(wp, 'mu') && wp.node >= 1 && wp.node <= obj.N+1
							mu_wp = wp.mu(:);
							% Find non-NaN components to constrain
							idx_constrained = ~isnan(mu_wp);
							if any(idx_constrained)
								constraints = [constraints; [vars.mu(idx_constrained, wp.node) == mu_wp(idx_constrained)]:sprintf('Waypoint Mean (node %d)', wp.node)];
							end
						end
					end
				end
			end

		end

		function constraints = convex_ineq(obj, vars)
			constraints = [
				[norm( chol(obj.P_f, 'lower') \ vars.S(:,:,obj.N+1), 2) - 1 <= 0]:'Terminal Covariance'
				% cone(vec(chol(obj.P_f, 'lower') \ vars.S(:,:,obj.N+1)), 1)
			];

			% Ensure positive diagonal elements of S for uniqueness
			for k = 1:obj.N+1
				constraints = [constraints
					[diag(vars.S(:,:,k)) >= 0]:'Positive Diagonal of S'
				];
			end

			% Add state chance constraints (if provided)
			if ~isempty(obj.chance_constraints_state)
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
					
					if strcmp(cc.type, 'affine')
						% Affine chance constraint: P(alpha'*x <= beta) >= 1-p
						% Note: affine chance constraints require mean variables
						if ~isfield(vars, 'mu') || isempty(vars.mu)
							error('State affine chance constraints require mean variables. Provide mu_0, mu_f, waypoints, or ensure mean variables are defined.');
						end
						alpha = cc.alpha(:); % ensure column vector
						beta = cc.beta;
						p = cc.p;
						z = norminv(1 - p);
						
						for k = nodes
							if k >= 1 && k <= obj.N+1
								constraints = [constraints;
									% alpha' * vars.mu(:,k) + z * norm(alpha' * vars.S(:,:,k)) - beta <= 0
									cone([ (- alpha' * vars.mu(:,k) + beta) / z; vars.S(:,:,k)' * alpha])
								];
							end
						end
					elseif strcmp(cc.type, 'norm')
						% Norm chance constraint: P(||x||_2 <= gamma) >= 1-p
						% Note: norm chance constraints require mean variables
						if ~isfield(vars, 'mu') || isempty(vars.mu)
							error('State norm chance constraints require mean variables. Provide mu_0, mu_f, waypoints, or ensure mean variables are defined.');
						end
						gamma = cc.gamma;
						p = cc.p;
						n = cc.n;
						q = sqrt(chi2inv(1 - p, n));
						
						for k = nodes
							if k >= 1 && k <= obj.N+1
								mu_k = vars.mu(:,k);
								S_k = vars.S(:,:,k);
								constraints = [constraints;
									[norm(mu_k, 2) + q * norm(S_k, 2) - gamma <= 0]:sprintf('State Norm Chance Constraint (k=%d)', k)
								];
							end
						end
					end
				end
			end

			% Add control chance constraints (if provided)
			if ~isempty(obj.chance_constraints_control)
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
					
					if strcmp(cc.type, 'affine')
						% Affine chance constraint: P(alpha'*u <= beta) >= 1-p
						% Note: affine chance constraints require mean variables
						if ~isfield(vars, 'v') || isempty(vars.v)
							error('Control affine chance constraints require mean variables. Provide mu_0, mu_f, waypoints, or ensure mean variables are defined.');
						end
						alpha = cc.alpha(:); % ensure column vector
						beta = cc.beta;
						p = cc.p;
						z = norminv(1 - p);
						
						for k = nodes
							v_k = vars.v(:,k);
							L_k = vars.L(:,:,k);
							constraints = [constraints;
								%[alpha' * v_k + z * norm(L_k'*alpha) - beta <= 0]:sprintf('Control Affine Chance Constraint (k=%d)', k)
								cone([ (- alpha' * v_k + beta) / z; L_k'*alpha])
								];
						end
					elseif strcmp(cc.type, 'norm')
						% Norm chance constraint: P(||u||_2 <= gamma) >= 1-p
						% Note: norm chance constraints require mean variables
						if ~isfield(vars, 'v') || isempty(vars.v)
							error('Control norm chance constraints require mean variables. Provide mu_0, mu_f, waypoints, or ensure mean variables are defined.');
						end
						gamma = cc.gamma;
						p = cc.p;
						n = cc.n;
						chi2q = chi2inv(1 - p, n);
						q = sqrt(chi2q);
						
						for k = nodes
							if k >= 1 && k <= obj.N
								v_k = vars.v(:,k);
								L_k = vars.L(:,:,k);
								constraints = [constraints;
									[norm(v_k, 2) + q * norm(L_k, 2) - gamma <= 0]:sprintf('Control Norm Chance Constraint (k=%d)', k)
								];
							end
						end
					end
				end

			end

		end

		function constraintLHS = noncvx_eq(obj, vars)
			% Non-convex equality constraints for QR-based covariance dynamics
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
				
				X_k = [A_k * S_k + B_k * L_k, G_k];

				[~, R_kp1_predicted] = obj.economy_qr_with_positive_diagonal(X_k');

				constraintLHS = [constraintLHS
					obj.vec_tril(S_kp1 - R_kp1_predicted')
				];
			end
		end

		function constraintLHS = noncvx_eq_relaxed(obj, vars, ref_vars)
			% Relaxed non-convex equality constraints using QR derivative linearization
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
				
				[dR, ~, R_X_k_ref] = obj.d_QR(X_k_ref', dX');

				constraintLHS = [constraintLHS
					obj.vec_tril(S_kp1 - R_X_k_ref' - dR')
				];

			end
		end

		% Non-convex inequality constraints
		function constraintLHS = noncvx_ineq(obj, vars)
			constraintLHS = [];
		end

		% Non-convex inequality constraints (relaxed)
		function constraintLHS = noncvx_ineq_relaxed(obj, vars, ref_vars)
			constraintLHS = [];
		end

		function constraints = convexified_exact(obj, vars, ref_vars)
			% Convex/convexified constraints that are imposed exactly but change with the reference variables
			% WARNING: In general, this breaks the convergence guarantee of SCvx/SCvx* and causes Delta L to become negative
			% Use this for handling circular keep-out zones which are linearized around the reference mean and covariance
			constraints = [];

			pos_idx = 1:2;

			for i = 1:length(obj.circular_obstacles)
				obstacle = obj.circular_obstacles{i};
				center = obstacle.center;
				radius = obstacle.radius;
				p = obstacle.p;
				z = norminv(1 - p);
				for k = 1:obj.N+1
					mu_k = vars.mu(:,k);
					S_k = vars.S(:,:,k);
					mu_ref_k = ref_vars.mu(:,k);
					a = - (mu_ref_k(pos_idx) - center);
					b = - 0.5 * norm(mu_ref_k(pos_idx) - center)^2  + 0.5 * radius^2 - a' * mu_ref_k(pos_idx);
					constraints = [constraints;
						cone([ (- a' * mu_k(pos_idx) - b) / z; S_k(pos_idx,pos_idx)' * a])
						% a' * mu_k(pos_idx) + z * norm(a' * S_k(pos_idx,pos_idx)) + b <= 0
					];
				end
			end
		end

		function postprocess(obj)
			obj.P = zeros(obj.nx, obj.nx, obj.N+1);
			obj.P_u = zeros(obj.nu, obj.nu, obj.N);
			obj.K = zeros(obj.nu, obj.nx, obj.N);
			for k = 1:obj.N+1
				S_k = value(obj.sol.S(:,:,k));
				obj.P(:,:,k) = S_k * S_k';
			end
			for k = 1:obj.N
				L_k = value(obj.sol.L(:,:,k));
				S_k = value(obj.sol.S(:,:,k));
				obj.P_u(:,:,k) = L_k * L_k';
				obj.K(:,:,k) = L_k / S_k;
			end
			% Mean variables (only if they were defined)
			if isfield(obj.sol, 'mu') && ~isempty(obj.sol.mu)
				obj.mu = value(obj.sol.mu);
			else
				obj.mu = zeros(obj.nx, obj.N+1);
			end
			if isfield(obj.sol, 'v') && ~isempty(obj.sol.v)
				obj.v = value(obj.sol.v);
			else
				obj.v = zeros(obj.nu, obj.N);
			end
		end

		function [dR, Qx, Rx] = d_QR(obj, X, dX, Qx, Rx, Rx_inv)
			arguments
				obj
				X
				dX
				Qx = []
				Rx = []
				Rx_inv = []
			end
			%D_QR Computes the differentials of the QR decomposition.
			% The dX input can be a sdpvar
			% If Qx and Rx are not provided, they are computed

			% Perform economy-size QR decomposition of the original matrix X
			if isempty(Qx) || isempty(Rx)
				[Qx, Rx] = obj.economy_qr_with_positive_diagonal(X);
			end

			if isempty(Rx_inv)
				Rx_inv = inv(Rx);
			end

			% Calculate the intermediate matrix V
			V = (Qx' * dX) * Rx_inv;

			% Create the anti-symmetric matrix 'A' from the lower triangular part of V.
			M = tril(V);
			A = M - M';

			dR = (V - A) * Rx;

		end

		function loss = compute_covariance_propagation_loss(obj)
			% Compute loss in covariance propagation according to:
			% Loss_k = ||φ(Kk, Pk) - Pk+1||F / ||Pk+1||F
			% where φ(Kk, Pk) := (Ak + Bk*Kk)*Pk*(Ak + Bk*Kk)' + Gk*Gk'
			%
			% Returns:
			%   loss: (N x 1) vector of losses for each time step k = 1, ..., N
			
			if isempty(obj.P) || isempty(obj.K)
				error('SqrtQRCovarianceSteering: P and K must be set. Solve the problem first.');
			end
			
			N = obj.N;
			loss = zeros(N, 1);
			
			for k = 1:N
				% Get system matrices for time step k
				if size(obj.A_sys, 3) == 1
					A_k = obj.A_sys;
				else
					A_k = obj.A_sys(:,:,k);
				end
				
				if size(obj.B_sys, 3) == 1
					B_k = obj.B_sys;
				else
					B_k = obj.B_sys(:,:,k);
				end
				
				if size(obj.G_sys, 3) == 1
					G_k = obj.G_sys;
				else
					G_k = obj.G_sys(:,:,k);
				end
				
				% Get feedback gain and covariance
				K_k = obj.K(:,:,k);
				P_k = obj.P(:,:,k);
				P_kp1 = obj.P(:,:,k+1);
				
				% Compute φ(Kk, Pk) = (Ak + Bk*Kk)*Pk*(Ak + Bk*Kk)' + Gk*Gk'
				A_closed = A_k + B_k * K_k;
				phi_Kk_Pk = A_closed * P_k * A_closed' + G_k * G_k';
				
				% Compute loss: ||φ(Kk, Pk) - Pk+1||F / ||Pk+1||F
				diff = phi_Kk_Pk - P_kp1;
				loss(k) = norm(diff, 'fro') / norm(P_kp1, 'fro');
			end
		end
				
	end

	methods (Static)

		function [Q, R] = economy_qr_with_positive_diagonal(X, options)
			arguments
				X
				options.check_rank = true
			end
			if options.check_rank
				if rank(X) < min(size(X))
					error('Input matrix X is rank deficient for QR decomposition.');
				end
			end

			[Q, R] = qr(X, "econ");

			% For uniqueness and consistent comparison, enforce the convention that R
			% has positive diagonal elements. This fixes sign ambiguities in Q and R.
			signs = diag(sign(diag(R)));
			Q = Q * signs;
			R = signs * R;
		end

		function out = vec_tril(M)
			% Returns the vectorized lower triangular part of matrix M
			n = size(M,1);
			idx = tril(true(n));
			out = M(idx);
		end
	end
end
