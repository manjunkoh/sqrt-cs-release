% Square root covariance steering optimizer with YALMIP optimizer capability
% This class inherits from SqrtQRCovarianceSteering and uses YALMIP's optimizer
% functionality to accelerate the algorithm by pre-compiling the optimization problem
classdef SqrtQRCovarianceSteeringOptimizer < SqrtQRCovarianceSteering

	properties 
		changing_params = struct();
		is_used_in_constraint_update = struct('S', true, 'L', true, 'mu', false, 'v', false);
	end

	methods
		function obj = SqrtQRCovarianceSteeringOptimizer(init_guess_struct, varargin)
			% Call parent constructor with all arguments
			obj@SqrtQRCovarianceSteering(init_guess_struct, varargin{:});
		end

		function val = get_use_optimizer(obj)
			val = true;
		end

		function set_changing_parameters_sdp(obj)
			% Define SDP variables for parameters that change between iterations
			% These correspond to the linearization points for the QR decomposition
			
			% Parameters for QR decomposition linearization
			obj.changing_params.S_ref = sdpvar(obj.nx, obj.nx, obj.N+1, 'full');
			obj.changing_params.L_ref = sdpvar(obj.nu, obj.nx, obj.N, 'full');
			
			% For each time step, we need the reference matrices for linearization
			for k = 1:obj.N+1
				obj.changing_params.S_ref(:,:,k) = tril(obj.changing_params.S_ref(:,:,k));
			end
		end

		function p = get_changing_parameters(obj, ref_vars)
			% Numerical update of the parameters that appear in the constraints
			% This function computes the reference values for linearization
			
			p.S_ref = ref_vars.S;
			p.L_ref = ref_vars.L;
		end

		function constraints = get_changing_constraints(obj, vars, ref_vars)
			% Get the changing constraints using YALMIP optimizer
			obj.set_changing_parameters_sdp();
			p = obj.changing_params;

			% Initialize constraints
			constraints = [];

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
				S_k_ref = p.S_ref(:,:,k);
				L_k_ref = p.L_ref(:,:,k);
				
				% Build the matrix for QR decomposition at reference point
				X_k_ref = [A_k * S_k_ref + B_k * L_k_ref, G_k];
				dX = [A_k * (S_k - S_k_ref) + B_k * (L_k - L_k_ref), zeros(size(G_k))];
				
				% Compute QR derivative using the d_QR function
				[dR, ~, R_X_k_ref] = d_QR(X_k_ref', dX');

				% Add linearized constraint (same structure as noncvx_eq_relaxed)
				constraintLHS = [constraintLHS
					obj.vec_tril(S_kp1 - R_X_k_ref' - dR')
				];
			end
			
			% The constraint is that the linearized constraint equals the slack variable
			constraints = constraintLHS == obj.slack_noncvx_eq;
		end

	end

end
