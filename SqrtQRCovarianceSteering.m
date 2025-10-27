% Square root covariance steering with QR decomposition-based covariance propagation
% This class inherits from SCPProblem 
classdef SqrtQRCovarianceSteering < SCPProblem

	properties
		init_guess_struct
		% imposing trust region on L is not recommended, as it causes chattering
		% near the solution and slower convergence.
		impose_trust_region_struct = struct('S', true, 'L', false, 'mu', false, 'v', false);
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
		chance_constraints_state = {}   % cell array (1 x N+1) of arrays of AffineChanceConstraint/NormChanceConstraint
		chance_constraints_control = {} % cell array (1 x N) of arrays of AffineChanceConstraint/NormChanceConstraint
		mu_0 = []  % initial mean (nx x 1) or empty -> assumed zero
		mu_f = []  % terminal mean (nx x 1) or empty -> assumed zero
		mu % state mean trajectory, set after solving
		v  % control mean trajectory, set after solving
	end

	methods
		function obj = SqrtQRCovarianceSteering(init_guess_struct, options)
			arguments
				init_guess_struct struct
				options.N
				options.nx
				options.nu
				options.nw
				options.A_sys
				options.B_sys
				options.G_sys
				options.P_0
				options.P_f
				options.Q
				options.R
				options.objective_type = 'LQG';
				options.chance_constraints_state = {}
				options.chance_constraints_control = {}
				options.mu_0 = []
				options.mu_f = []
			end
			obj@SCPProblem();
			obj.init_guess_struct = init_guess_struct;

			obj.N = options.N;
			obj.nx = options.nx;
			obj.nu = options.nu;
			obj.nw = options.nw;
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

			obj.initialize();
		end

		function vars = define_vars(obj)
			vars.S = sdpvar(obj.nx, obj.nx, obj.N+1);
			for k = 1:obj.N+1
				vars.S(:,:,k) = tril(vars.S(:,:,k));
			end
			vars.L = sdpvar(obj.nu, obj.nx, obj.N);
			vars.mu = sdpvar(obj.nx, obj.N+1);
			vars.v = sdpvar(obj.nu, obj.N);
		end

		function J0 = objective(obj, vars)
			% Minimize control effort and state covariance
			J0 = 0;
			
			switch obj.objective_type
				case 'LQG'
					% Add control cost
					for k = 1:obj.N
						J0 = J0 + trace(vars.L(:,:,k) * vars.L(:,:,k)' * obj.R);
					end
					
					% Add state covariance cost
					for k = 1:obj.N
						J0 = J0 + trace(vars.S(:,:,k) * vars.S(:,:,k)' * obj.Q);
					end

					% Add mean state and control cost
					for k = 1:obj.N
						J0 = J0 + vars.mu(:,k)' * obj.Q * vars.mu(:,k) + vars.v(:,k)' * obj.R * vars.v(:,k);
					end
				case 'DV99'
					% Add control cost only
					for k = 1:obj.N
						J0 = J0 + norm(vars.L(:,:,k), 2);
					end
			end
		end

		function constraints = convex_eq(obj, vars)
			constraints = [
				[vars.S(:,:,1) == chol(obj.P_0, 'lower')]:'Initial Covariance'
			];

			% initial mean if provided, otherwise zero
			if ~isempty(obj.mu_0)
				constraints = [constraints; [vars.mu(:,1) == obj.mu_0]:'Initial Mean'];
			end
			for k = 1:obj.N
				A_k = obj.A_sys(:,:,k);
				B_k = obj.B_sys(:,:,k);
				constraints = [constraints; [vars.mu(:,k+1) == A_k * vars.mu(:,k) + B_k * vars.v(:,k)]:'Mean Dynamics'];
			end
			% enforce terminal mean if provided in mu_f
			if ~isempty(obj.mu_f)
				constraints = [constraints; [vars.mu(:,obj.N+1) == obj.mu_f]:'Terminal Mean'];
			end
		end

		function constraints = convex_ineq(obj, vars)
			constraints = [
				[norm( chol(obj.P_f, 'lower') \ vars.S(:,:,obj.N+1), 2) - 1 <= 0]:'Terminal Covariance'
			];

			for k = 1:obj.N+1
				constraints = [constraints
					[diag(vars.S(:,:,k)) >= 0]:'Positive Diagonal of S'
				];
			end

			% for j = 1:length(obj.chance_constraints_state)
			% 	c = obj.chance_constraints_state(j);
			% 	constraints = [constraints; c.toYALMIPConstraint(vars.mu, vars.S)];
			% end
			
			% Add state chance constraints (if provided)
			p = 0.005;
			a = [0.2 -1 0 0]';
			b = 0.2;
			for k = 1:obj.N
				mu_k = vars.mu(:,k);
				S_k = vars.S(:,:,k);

				constraints = [constraints
					[norminv(1-p) * norm(a' * S_k) + a' * mu_k - b <= 0]:'State Chance Constraint'
					];
			end

			a = [0.2 1 0 0]';
			b = 0.2;
			for k = 1:obj.N
				mu_k = vars.mu(:,k);
				S_k = vars.S(:,:,k);

				constraints = [constraints
					[norminv(1-p) * norm(a' * S_k) + a' * mu_k - b <= 0]:'State Chance Constraint'
					];
			end

			% Add control chance constraints (if provided)
			% for k = 1:obj.N
			% 	v_k = vars.v(:,k);
			% 	L_k = vars.L(:,:,k);
			% 	for c = obj.chance_constraints_control
			% 		constraints = [constraints; c{1}.toYALMIPConstraint(v_k, L_k)];
			% 	end
			% end
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

				S_kp1_predicted = qr(X_k', "econ");
				% make R have positive diagonal
				signs = diag(sign(diag(S_kp1_predicted)));
				S_kp1_predicted = signs * S_kp1_predicted;

				S_kp1_predicted = S_kp1_predicted';

				constraintLHS = [constraintLHS
					obj.vec_tril(S_kp1 - S_kp1_predicted)
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
			obj.mu = value(obj.sol.mu);
			obj.v = value(obj.sol.v);
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
				
	end

	methods (Static)

		function [Q, R] = economy_qr_with_positive_diagonal(X)
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
