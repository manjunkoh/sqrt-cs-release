classdef CovarianceSteering < handle
	% Covariance Steering class
	% 
	% This class implements covariance steering optimization problems with two solution methods:
	% 1. 'full_covariance' - Original method using full covariance matrices (default)
	% 2. 'block_cholesky' - Block Cholesky method based on Okamoto et al. (2018)
	%
	% Usage:
	%   cs = CovarianceSteering('A', A, 'B', B, 'D', D, 'P_0', P_0, 'P_f', P_f, ...
	%                          'Q', Q, 'R', R, 'nx', nx, 'nu', nu, 'N', N, ...
	%                          'solution_method', 'block_cholesky');
	%   diagnostic = cs.solve_problem();
	properties
		A % System dynamics matrix
		B % Input matrix
		D % Process noise input matrix (Cholesky factor of process noise covariance)
		P_0 % Initial state covariance
		P_f % Final state covariance
		Q % State cost matrix; assume constant over time
		R % Input cost matrix; assume constant over time
		nx % State dimension
		nu % Input dimension
		nw % Process noise dimension
		N % Time horizon
		rho % Relaxation parameter for final covariance constraint
		std_cdn = 0; % standard deviation of control-dependent noise
		sdpvars % Struct to hold YALMIP sdpvar variables
		covariance_scaling = 1; % Scaling factor for state covariance
		K % Feedback gains; computed after solving
			% For 'full_covariance': K(:,:,k) is the feedback gain at time k
			% For 'block_cholesky': K(:,:,k,i) is the gain at time k for state at time i
		yalmip_settings = sdpsettings('verbose', 0, 'solver', 'mosek', 'savesolveroutput', true);
		u_mean = [];
		solution_method = 'full_covariance'; % Solution method: 'full_covariance' or 'block_cholesky'
		optimal_objective = NaN; % Optimal objective value after solving
	end
	
	methods
		function obj = CovarianceSteering(options)
			arguments
				options.A
				options.B
				options.D
				options.P_0
				options.P_f
				options.Q
				options.R
				options.nx
				options.nu
				options.N
				options.covariance_scaling = 1
				options.rho = []
				options.std_cdn = 0
				options.u_mean = []
				options.solution_method = 'full_covariance'
			end
			
			obj.A = options.A;
			obj.B = options.B;
			obj.D = options.D;
			obj.P_0 = options.P_0;
			obj.P_f = options.P_f;
			obj.Q = options.Q;
			obj.R = options.R;
			obj.nx = options.nx;
			obj.nu = options.nu;
			obj.nw = size(options.D, 2);
			obj.N = options.N;
			obj.covariance_scaling = options.covariance_scaling;
			obj.rho = options.rho;
			obj.std_cdn = options.std_cdn;
			obj.solution_method = options.solution_method;
			obj.sdpvars = struct('P', [], 'U', [], 'Y', [], 'J', [], 'constraints', []);
			if isempty(options.u_mean)
				obj.u_mean = zeros(obj.nu, obj.N);
			else
				assert(size(options.u_mean, 1) == obj.nu)
				assert(size(options.u_mean, 2) >= obj.N)
				obj.u_mean = options.u_mean;
			end

		end

		function set_sdpvars(obj)
			% Define YALMIP sdpvar variables
			nx = obj.nx;
			nu = obj.nu;
			N = obj.N;

			obj.sdpvars.P = sdpvar(nx, nx, N+1); % State covariance
			obj.sdpvars.U = sdpvar(nu, nx, N); % Control gain
			obj.sdpvars.Y = sdpvar(nu, nu, N); % Auxiliary variable for input covariance
		end

		function set_objective(obj)
			J = 0;
			for k = 1:obj.N
				J = J + trace(obj.Q * obj.sdpvars.P(:,:,k)) ...
					+ trace(obj.R * obj.sdpvars.Y(:,:,k));
			end
			obj.sdpvars.J = J;
		end

		function set_constraints(obj)
			% Set all constraints at once
			constraints = [
				obj.get_boundary_constraints()
				obj.get_dynamics_constraints()
				obj.get_LMI_constraints(1:obj.N)
				obj.get_positive_semidefinite_constraints()
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

		function constraints = get_dynamics_constraints(obj)
			A = obj.A;
			B = obj.B;
			D = obj.D;

			constraints = [];
			for k = 1:obj.N
				constraints = [constraints
					obj.sdpvars.P(:,:,k+1) == ...
						(A(:,:,k) * obj.sdpvars.P(:,:,k) * A(:,:,k)' ...
						+ A(:,:,k) * obj.sdpvars.U(:,:,k)' * B(:,:,k)' ...
						+ B(:,:,k) * obj.sdpvars.U(:,:,k) * A(:,:,k)' ...
						+ (1 + obj.std_cdn^2) * B(:,:,k) * obj.sdpvars.Y(:,:,k) * B(:,:,k)' ...
						+ obj.std_cdn^2 * B(:,:,k) * obj.u_mean(:,k) * obj.u_mean(:,k)' * B(:,:,k)' ...
						+ (D(:,:,k) * D(:,:,k)') .* obj.covariance_scaling)
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
					% [obj.sdpvars.P(:,:,k) >= 0]:'P_positive_semidefinite'
				];
			end
		end

		function diagnostic = solve_problem(obj, options)
			% Solve the SDP using YALMIP
			arguments
				obj
				options.verbose = 0
				options.solver = 'mosek'
				options.savesolveroutput = true
			end
			
			% fprintf('Solving covariance steering problem with N=%d using %s method...\n', obj.N, obj.solution_method);
			
			% Choose solution method
			if strcmp(obj.solution_method, 'full_covariance')
				obj.set_sdpvars();
				obj.set_objective();
				obj.set_constraints();
			elseif strcmp(obj.solution_method, 'block_cholesky')
				obj.set_sdpvars_block_cholesky();
				obj.set_objective_block_cholesky();
				obj.set_constraints_block_cholesky();
			else
				error('Unknown solution method: %s. Use ''full_covariance'' or ''block_cholesky''', obj.solution_method);
			end

			% obj.sdpvars.constraints
			settings = obj.yalmip_settings;
			settings.verbose = options.verbose;
			settings.solver = options.solver;
			settings.savesolveroutput = options.savesolveroutput;

			diagnostic = optimize(obj.sdpvars.constraints, obj.sdpvars.J, settings);

			switch diagnostic.problem
				case 0
					fprintf('Optimization successful. Objective value: %f\n', value(obj.sdpvars.J));
					if strcmp(obj.solution_method, 'full_covariance')
						obj.set_feedback_gains();
						obj.check_lossless(value(obj.sdpvars.P), obj.K, obj.A, obj.B, obj.D, obj.N);
					elseif strcmp(obj.solution_method, 'block_cholesky')
						obj.set_feedback_gains_block_cholesky();
					end
					obj.optimal_objective = value(obj.sdpvars.J);
				case 1
					warning('YALMIP:Infeasible', 'Problem infeasible');
				case 4
					if strcmp(obj.solution_method, 'full_covariance')
						obj.set_feedback_gains();
						fprintf("Numerical problems. Objective value: %f.\n", value(obj.sdpvars.J));
						obj.check_lossless(value(obj.sdpvars.P), obj.K, obj.A, obj.B, obj.D, obj.N);
					elseif strcmp(obj.solution_method, 'block_cholesky')
						obj.set_feedback_gains_block_cholesky();
						% warning("Numerical problems. Objective value: %f. Check https://themosekblog.blogspot.com/2014/06/what-if-solver-stall.html for details.", value(obj.sdpvars.J));
					end
					obj.optimal_objective = value(obj.sdpvars.J);
				otherwise
					warning('YALMIP:Failed', 'Solver failed with code %d', diagnostic.problem);
			end
		end

		function set_feedback_gains(obj)
			% Compute feedback gains from optimal U and P
			nx = obj.nx;
			nu = obj.nu;
			N = obj.N;
			K = zeros(nu, nx, N);
			for k = 1:N
				P_k = value(obj.sdpvars.P(:,:,k));
				U_k = value(obj.sdpvars.U(:,:,k));
				K(:,:,k) = U_k / P_k; % Feedback gain
			end
			obj.K = K;
        end

        function [is_lossless, worst_loss] = check_lossless(obj, P, K, A, B, D, N, options)
			arguments
				obj
				P
				K
				A
				B
				D
				N
				options.tol = 1E-4
			end
				
			worst_loss = 0;
			is_lossless = true;
            for k = 1:N
                P_next = (A(:,:,k) + B(:,:,k)*K(:,:,k)) * P(:,:,k) * (A(:,:,k) + B(:,:,k)*K(:,:,k))' + D(:,:,k)*D(:,:,k)';
				% Add control-dependent noise term
				P_next = P_next + obj.std_cdn^2 * B(:,:,k) * K(:,:,k) * P(:,:,k) * K(:,:,k)' * B(:,:,k)';
				P_next = P_next + obj.std_cdn^2 * B(:,:,k) * obj.u_mean(:,k) * obj.u_mean(:,k)' * B(:,:,k)';
                % assert(norm(P_next - P(:,:,k+1), 'fro') < options.tol)
				loss = norm(P_next - P(:,:,k+1), 'fro') / norm(P(:,:,k+1), 'fro');
				if loss > options.tol
					is_lossless = false;
					worst_loss = max(worst_loss, loss);
				end
            end

			if is_lossless
				fprintf('Losslessness verified within tolerance %g.\n', options.tol);
			else
				fprintf('Losslessness NOT verified. Worst relative loss\n: %g', worst_loss);
			end
        end

		function varargout = plot_Y_lambda_max(obj, varargin)
			N = obj.N;
			
			P_u = obj.get_control_covariance();
			
			% Compute maximum eigenvalues of control covariance matrices
			Y_lmd_max = NaN(N, 1);
			for k = 1:N
				% Ensure matrix is symmetric for eigenvalue computation
				P_u_k = (P_u(:,:,k) + P_u(:,:,k)') / 2;
				Y_lmd_max(k) = lambda_max(P_u_k);
			end
			
			[varargout{1:nargout}] = stairsZOH(0:N, Y_lmd_max, varargin{:});
			xlabel('Node $k$')
			ylabel('$\lambda_{\max}(Y_k)$')
		end

		function plot_state_sigmas(obj, varargin)
			nx = obj.nx;
			N = obj.N;
			
			P = obj.get_state_covariance();
			
			% Compute standard deviations from covariance matrices
			sigmas = NaN(nx, N+1);
			for k = 1:N+1
				for i = 1:nx
					sigmas(i,k) = sqrt(P(i,i,k));
				end
			end

			hold on
			colors = lines(nx);
			for i = 1:nx
				plot(0:N, sigmas(i,:), 'Color', colors(i,:), 'DisplayName', sprintf('$\\sigma_{x_%d}$', i), varargin{:});
			end
			hold off
			xlabel('Node $k$')
			ylabel('State standard deviations')
			legend('Location', 'best')
		end

		% Block Cholesky solution method
		function set_sdpvars_block_cholesky(obj)
			% Define YALMIP sdpvar variables for block Cholesky method
			nx = obj.nx;
			nu = obj.nu;
			N = obj.N;

			% F is the block diagonal feedback gain matrix
			% obj.sdpvars.F = [obj.block_diag_sdp(nu, nx, N) zeros(N*nu, nx)];
			obj.sdpvars.F = [obj.block_tril_sdp(nu, nx, N) zeros(N*nu, nx)];
		end

		function set_objective_block_cholesky(obj)
			% Set objective for block Cholesky method
			[A_blk, B_blk, ~, ~, G_blk] = obj.get_block_matrices();
			
			% Create block diagonal Q and R matrices
			Q_blk = obj.repblkdiag(obj.Q, obj.N+1);
			R_blk = obj.repblkdiag(obj.R, obj.N);
			
			% Compute S matrix
			S = A_blk * obj.P_0 * A_blk' + G_blk * G_blk';
			I = eye((obj.N+1) * obj.nx);
			
			% Objective function from Okamoto et al (2018)
			obj.sdpvars.J = trace(((I + B_blk * obj.sdpvars.F)' * Q_blk * (I + B_blk * obj.sdpvars.F) + obj.sdpvars.F' * R_blk * obj.sdpvars.F) * S);
		end

		function set_constraints_block_cholesky(obj)
			% Set constraints for block Cholesky method
			[~, B_blk, Ex, ~, ~, chol_S] = obj.get_block_matrices();
			
			I = eye((obj.N+1) * obj.nx);
			
			% Terminal constraint
			constraints = [
				norm(chol(obj.P_f, 'lower') \ (Ex(:,:,obj.N+1) * (I + B_blk * obj.sdpvars.F) * chol_S)) - 1 <= 0
				% norm(reshape(chol(obj.P_f, 'lower') \ (Ex(:,:,obj.N+1) * (I + B_blk * obj.sdpvars.F) * chol_S), [], 1)) - 1 <= 0
			];
			
			obj.sdpvars.constraints = constraints;
		end

		function set_feedback_gains_block_cholesky(obj)
			% Compute feedback gains from block Cholesky solution
			[~, B_blk, ~, ~, ~, ~] = obj.get_block_matrices();
			I = eye((obj.N+1) * obj.nx);
			
			% Extract block diagonal feedback gains
			F = value(obj.sdpvars.F);
			K_blk = F / (I + B_blk * F);
			
			% Convert to time-indexed gains
			nx = obj.nx;
			nu = obj.nu;
			N = obj.N;
			K = zeros(nu, nx, N, N);
			
			for k = 1:N
				for i = 1:k
					K(:,:,k,i) = K_blk(1+(k-1)*nu:k*nu, 1+(i-1)*nx:i*nx);
				end
			end
			
			obj.K = K;
		end

		% Block Cholesky utility functions
		function [A_blk, B_blk, Ex, Eu, G_blk, chol_S] = get_block_matrices(obj)
			% Construct block matrices for the block Cholesky method
			nx = obj.nx;
			nu = obj.nu;
			nw = obj.nw;
			N = obj.N;
			A_sys = obj.A;
			B_sys = obj.B;
			G_sys = obj.D; % D is the process noise input matrix
			P_0 = obj.P_0;

			if size(A_sys, 3) == 1
				A_sys = repmat(A_sys, [1 1 N]);
				B_sys = repmat(B_sys, [1 1 N]);
				G_sys = repmat(G_sys, [1 1 N]);
			end

			A_blk = zeros((N+1)*nx, nx);
			A_blk(1:nx,:) = eye(nx);
			B_blk = zeros((N+1)*nx, N*nu);
			G_blk = zeros((N+1)*nx, N*nw);
			B_blk(nx+1:2*nx, 1:nu) = B_sys(:,:,1);
			G_blk(nx+1:2*nx, 1:nw) = G_sys(:,:,1);

			for i = 1:N
				A_blk(i*nx+1:(i+1)*nx, :) = A_sys(:,:,i) * A_blk((i-1)*nx+1:i*nx, :);
			end
			for i = 2:N
				B_blk(i*nx+1:(i+1)*nx, 1:i*nu) = [A_sys(:,:,i) * B_blk((i-1)*nx+1:i*nx, 1:(i-1)*nu) B_sys(:,:,i)];
				G_blk(i*nx+1:(i+1)*nx, 1:i*nw) = [A_sys(:,:,i) * G_blk((i-1)*nx+1:i*nx, 1:(i-1)*nw) G_sys(:,:,i)];
			end

			Eu = zeros(nu, N*nu, N);
			Ex = zeros(nx, (N+1)*nx, N+1);
			for k = 1:N
				Eu(:,:,k) = [zeros(nu,(k-1)*nu)  eye(nu) zeros(nu,(N-k)*nu)];
			end
			for k = 1:N+1
				Ex(:,:,k) = [zeros(nx,(k-1)*nx) eye(nx) zeros(nx,(N-k+1)*nx)];
			end
			chol_S = [A_blk * chol(P_0, 'lower') G_blk];
		end

		function Q_blk = repblkdiag(obj, M, n)
			% Create block diagonal matrix by repeating M n times
			C = repmat({M}, n, 1);
			Q_blk = blkdiag(C{:});
		end

		function L = block_diag_sdp(obj, x, y, N)
			% Returns a (Nx,Ny) sdpvar of a (x,y) block N by N diagonal matrix
			% More efficient than naively indexing the matrix to set 0 constraints
			L_block = sdpvar(repmat(x,1,N), repmat(y,1,N));
			L = blkdiag(L_block{:});
		end

		function L = block_tril_sdp(obj, x,y,N)
			% returns a (Nx,Ny) sdpvar of a (x,y) block N by N lower triangular matrix
			% more efficient than naively indexing the matrix to set 0 constraints
			% https://groups.google.com/g/yalmip/c/PHnGnHXNcEA
			L = kron(tril(ones(N)),ones(x,y));
			L = L.*sdpvar(N*x,N*y);
		end

		function P = get_state_covariance(obj)
			% Get state covariance matrices for both solution methods
			nx = obj.nx;
			N = obj.N;
			
			switch obj.solution_method
				case 'full_covariance'
					P = value(obj.sdpvars.P);
				case 'block_cholesky'
					% For block Cholesky method, compute state covariances from block matrices
					[A_blk, B_blk, Ex, ~, G_blk, ~] = obj.get_block_matrices();
					F = value(obj.sdpvars.F);
					I = eye((obj.N+1) * obj.nx);
					
					% Compute the full state covariance matrix
					P_X = (I + B_blk * F) * (A_blk * obj.P_0 * A_blk' + G_blk * G_blk') * (I + B_blk * F)';
					
					% Extract state covariances at each time step
					P = zeros(nx, nx, N+1);
					for k = 1:N+1
						P(:,:,k) = Ex(:,:,k) * P_X * Ex(:,:,k)';
					end
			end
		end

		function P_u = get_control_covariance(obj)
			% Get control covariance matrices for both solution methods
			nu = obj.nu;
			N = obj.N;
			
			switch obj.solution_method
				case 'full_covariance'
					P_u = value(obj.sdpvars.Y);
				case 'block_cholesky'
					F = value(obj.sdpvars.F);
					[A_blk, B_blk, ~, Eu, G_blk, ~] = obj.get_block_matrices();
					
					% Compute S matrix and P_U (input covariance)
					S = A_blk * obj.P_0 * A_blk' + G_blk * G_blk';
					P_U = F * S * F';
					
					% Extract control covariances at each time step
					P_u = zeros(nu, nu, N);
					for k = 1:N
						P_u(:,:,k) = Eu(:,:,k) * P_U * Eu(:,:,k)';
					end
			end
		end

		function J_step = get_objective(obj)
			% Get objective function value per step for both solution methods
			nx = obj.nx;
			nu = obj.nu;
			N = obj.N;
			
			J_step = zeros(N, 1);
			
			P = obj.get_state_covariance();
			P_u = obj.get_control_covariance();
			
			for k = 1:N
				% State cost: trace(Q * P_k)
				state_cost = trace(obj.Q * P(:,:,k));
				% Control cost: trace(R * P_u_k)
				control_cost = trace(obj.R * P_u(:,:,k));
				J_step(k) = state_cost + control_cost;
			end
					
		end

	end
end