classdef BlockCholeskySteering < CovarianceSteeringBase
	% Block Cholesky Covariance Steering class
	% 
	% This class implements covariance steering using the block Cholesky method
	% based on Okamoto et al. (2019).
	%
	% Usage:
	%   cs = BlockCholeskySteering('A', A, 'B', B, 'G', G, 'P_0', P_0, 'P_f', P_f, ...
	%                              'Q', Q, 'R', R, 'nx', nx, 'nu', nu, 'N', N);
	%   diagnostic = cs.solve_problem();
	
	methods
		function obj = BlockCholeskySteering(varargin)
			obj = obj@CovarianceSteeringBase(varargin{:});
		end

		function set_sdpvars(obj)
			% Define YALMIP sdpvar variables for block Cholesky method
			% Following Okamoto et al. (2019): V is open-loop control mean, K is block diagonal
			nx = obj.nx;
			nu = obj.nu;
			N = obj.N;

			% K is the block diagonal feedback gain matrix (Okamoto 2019)
			obj.sdpvars.K = [obj.block_diag_sdp(nu, nx, N) zeros(N*nu, nx)];
			
			% V is the open-loop control mean (vectorized) - needed if chance constraints or mean boundary conditions
			if ~isempty(obj.chance_constraints_state) || ~isempty(obj.chance_constraints_control) ...
					|| ~isempty(obj.mu_0) || ~isempty(obj.mu_f)
				obj.sdpvars.V = sdpvar(N*nu, 1, 'full');
			else
				obj.sdpvars.V = [];
			end
			
			obj.sdpvars.J = [];
			obj.sdpvars.constraints = [];
		end

		function set_objective(obj)
			switch obj.objective_type
				case {'LQG', 'LQR', 'LQ'}
					% Set objective for block Cholesky method
					[A_blk, B_blk, ~, ~, G_blk, ~] = obj.get_block_matrices();
					
					% Create block diagonal Q and R matrices
					Q_blk = obj.repblkdiag(obj.Q, obj.N+1);
					R_blk = obj.repblkdiag(obj.R, obj.N);
					
					% Compute S matrix
					S = A_blk * obj.P_0 * A_blk' + G_blk * G_blk';
					I = eye((obj.N+1) * obj.nx);
					
					% Covariance part of objective (always included) - using K (Okamoto 2019)
					J_cov = trace(((I + B_blk * obj.sdpvars.K)' * Q_blk * (I + B_blk * obj.sdpvars.K) + obj.sdpvars.K' * R_blk * obj.sdpvars.K) * S);
					
					% Mean part of objective (if V is defined, following Okamoto et al. 2019)
					if ~isempty(obj.sdpvars.V)
						% Initialize x_0_bar
						if ~isempty(obj.mu_0)
							x_0_bar = obj.mu_0(:);
						else
							x_0_bar = zeros(obj.nx, 1);
						end
						
						% Mean part: (A_blk*x_0_bar+B_blk*V)'*Q_blk*(A_blk*x_0_bar+B_blk*V) + V'*R_blk*V
						J_mean = (A_blk * x_0_bar + B_blk * obj.sdpvars.V)' * Q_blk * (A_blk * x_0_bar + B_blk * obj.sdpvars.V) + obj.sdpvars.V' * R_blk * obj.sdpvars.V;
						
						obj.sdpvars.J = J_mean + J_cov;
					else
						obj.sdpvars.J = J_cov;
					end
					
				case 'DV99'
					% DV99 objective: minimize 99th percentile of control norm
					% Formulation: sum_k [norm(Eu_k * V) + q * norm(Eu_k * K * chol_S)]
					% where q = sqrt(chi2inv(0.99, nu))
					% Similar to norm control chance constraints
					if isempty(obj.sdpvars.V)
						error('BlockCholeskySteering: Mean control variables (V) must be defined for DV99 objective.');
					end
					
					[~, ~, ~, Eu, ~, chol_S] = obj.get_block_matrices();
					q = sqrt(chi2inv(0.99, obj.nu));
					obj.sdpvars.J = 0;
					
					for k = 1:obj.N
						% Mean part: norm(Eu_k * V)
						V_bar_k = Eu(:,:,k) * obj.sdpvars.V;
						
						% Covariance part: q * norm(Eu_k * K * chol_S)
						L_k = Eu(:,:,k) * obj.sdpvars.K * chol_S;
						
						obj.sdpvars.J = obj.sdpvars.J + norm(V_bar_k, 2) + q * norm(L_k, 2);
					end
					
				otherwise
					error('BlockCholeskySteering: Unknown objective type ''%s''.', obj.objective_type);
			end
		end

		function set_constraints(obj)
			% Set constraints for block Cholesky method
			% Following Okamoto et al. (2019) formulation
			[A_blk, B_blk, Ex, Eu, ~, chol_S] = obj.get_block_matrices();
			
			I = eye((obj.N+1) * obj.nx);
			
			% Initialize x_0_bar (initial mean)
			if ~isempty(obj.mu_0)
				x_0_bar = obj.mu_0(:);
			else
				x_0_bar = zeros(obj.nx, 1);
			end
			
			% Terminal covariance constraint - using K (Okamoto 2019)
			obj.sdpvars.constraints = [
				norm(chol(obj.P_f, 'lower') \ (Ex(:,:,obj.N+1) * (I + B_blk * obj.sdpvars.K) * chol_S)) - 1 <= 0
			];
			
			% Terminal mean constraint (if specified)
			if ~isempty(obj.mu_f) && ~isempty(obj.sdpvars.V)
				x_f_bar = obj.mu_f(:);
				obj.sdpvars.constraints = [
					obj.sdpvars.constraints;
					Ex(:,:,obj.N+1) * (A_blk * x_0_bar + B_blk * obj.sdpvars.V) == x_f_bar
				];
			end
			
			% State chance constraints (following Okamoto et al. 2019)
			if ~isempty(obj.chance_constraints_state) && ~isempty(obj.sdpvars.V)
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
						% Following Okamoto 2019: norminv(1-p) * norm(alpha' * E * (I+B*K) * chol_S) 
						%                        + alpha' * E * (A*x_0 + B*V) <= beta
						alpha = cc.alpha(:); % ensure column vector
						beta = cc.beta;
						p = cc.p;
						z = norminv(1 - p);
						
						for k = nodes
							if k >= 1 && k <= obj.N+1
								% Mean part: alpha' * E_k * (A_blk * x_0_bar + B_blk * V)
								mean_part = alpha' * Ex(:,:,k) * (A_blk * x_0_bar + B_blk * obj.sdpvars.V);
								
								% Covariance part: z * norm(alpha' * E_k * (I + B_blk * K) * chol_S)
								cov_part = z * norm(alpha' * Ex(:,:,k) * (I + B_blk * obj.sdpvars.K) * chol_S);
                                
								obj.sdpvars.constraints = [
									obj.sdpvars.constraints;
									cov_part + mean_part - beta <= 0
								];

							end
						end
					elseif strcmp(cc.type, 'norm')
						% Norm chance constraint: P(||x||_2 <= gamma) >= 1-p
						gamma = cc.gamma;
						p = cc.p;
						n = cc.n;
						q = sqrt(chi2inv(1 - p, n));
						
						for k = nodes
							if k >= 1 && k <= obj.N+1
								% Mean part: norm(E_k * (A_blk * x_0_bar + B_blk * V))
								X_bar_k = Ex(:,:,k) * (A_blk * x_0_bar + B_blk * obj.sdpvars.V);
								
								% Covariance part: q * norm(E_k * (I + B_blk * K) * chol_S)
								S_k = Ex(:,:,k) * (I + B_blk * obj.sdpvars.K) * chol_S;
								
								obj.sdpvars.constraints = [
									obj.sdpvars.constraints;
									norm(X_bar_k, 2) + q * norm(S_k, 2) - gamma <= 0
								];
							end
						end
					end
				end
			end
			
			% Control chance constraints (following Okamoto et al. 2019)
			if ~isempty(obj.chance_constraints_control) && ~isempty(obj.sdpvars.V)
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
						% Following Okamoto 2019: norminv(1-p) * norm(alpha' * Eu * K * chol_S) 
						%                        + alpha' * Eu * V <= beta
						alpha = cc.alpha(:); % ensure column vector
						beta = cc.beta;
						p = cc.p;
						z = norminv(1 - p);
						
						for k = nodes
							if k >= 1 && k <= obj.N
								% Mean part: alpha' * Eu_k * V
								mean_part = alpha' * Eu(:,:,k) * obj.sdpvars.V;
								
								% Covariance part: z * norm(alpha' * Eu_k * K * chol_S)
								cov_part = z * norm(alpha' * Eu(:,:,k) * obj.sdpvars.K * chol_S);
								
								obj.sdpvars.constraints = [
									obj.sdpvars.constraints;
									cov_part + mean_part - beta <= 0
								];
							end
						end
						
					elseif strcmp(cc.type, 'norm')
						% Norm chance constraint: P(||u||_2 <= gamma) >= 1-p
						gamma = cc.gamma;
						p = cc.p;
						n = cc.n;
						chi2q = chi2inv(1 - p, n);
						q = sqrt(chi2q);
						
						for k = nodes
							if k >= 1 && k <= obj.N
								% Mean part: norm(Eu_k * V)
								V_bar_k = Eu(:,:,k) * obj.sdpvars.V;
								
								% Covariance part: q * norm(Eu_k * K * chol_S)
								L_k = Eu(:,:,k) * obj.sdpvars.K * chol_S;
								
								obj.sdpvars.constraints = [
									obj.sdpvars.constraints;
									norm(V_bar_k, 2) + q * norm(L_k, 2) - gamma <= 0
								];
							end
						end
						
					else
						error('BlockCholeskySteering: Unsupported control chance constraint type.');
					end
				end
			end
		end

		function set_feedback_gains(obj)
			% Compute feedback gains from block Cholesky solution
			% Following Okamoto et al. (2019): K is block diagonal, extract diagonal blocks
			nx = obj.nx;
			nu = obj.nu;
			N = obj.N;
			
			% Extract block diagonal feedback gains (Okamoto 2019)
			K_blk = value(obj.sdpvars.K);
			
			% Convert to time-indexed gains (only diagonal blocks in Okamoto 2019)
			% K(:,:,k) is the feedback gain at time k
			K = zeros(nu, nx, N);
			
			for k = 1:N
				% Extract diagonal block: K_k = K_blk(k-th block)
				K(:,:,k) = K_blk(1+(k-1)*nu:k*nu, 1+(k-1)*nx:k*nx);
			end
			
			obj.K = K;
			
			% Compute and store state and control covariances
			[A_blk, B_blk, Ex, Eu, G_blk, ~] = obj.get_block_matrices();
			I = eye((obj.N+1) * obj.nx);
			
			% Compute S matrix
			S = A_blk * obj.P_0 * A_blk' + G_blk * G_blk';
			
			% Compute the full state covariance matrix - using K (Okamoto 2019)
			P_X = (I + B_blk * K_blk) * S * (I + B_blk * K_blk)';
			
			% Extract state covariances at each time step
			obj.P = zeros(nx, nx, N+1);
			for k = 1:N+1
				obj.P(:,:,k) = Ex(:,:,k) * P_X * Ex(:,:,k)';
			end
			
			% Compute control covariances
			P_U = K_blk * S * K_blk';
			obj.P_u = zeros(nu, nu, N);
			for k = 1:N
				obj.P_u(:,:,k) = Eu(:,:,k) * P_U * Eu(:,:,k)';
			end
			
			% Compute mean trajectories from V (following Okamoto et al. 2019)
			if ~isempty(obj.sdpvars.V)
				V = value(obj.sdpvars.V);
				
				% Initialize x_0_bar
				if ~isempty(obj.mu_0)
					x_0_bar = obj.mu_0(:);
				else
					x_0_bar = zeros(nx, 1);
				end
				
				% Compute state mean trajectory: X_bar = A_blk * x_0_bar + B_blk * V
				X_bar = A_blk * x_0_bar + B_blk * V;
				obj.mu = reshape(X_bar, nx, N+1);
				
				% Extract control mean trajectory: v = reshape(V, nu, N)
				obj.v = reshape(V, nu, N);
			else
				obj.mu = [];
				obj.v = [];
			end
		end

		% Block Cholesky utility functions
		function [A_blk, B_blk, Ex, Eu, G_blk, chol_S] = get_block_matrices(obj)
			% Construct block matrices for the block Cholesky method (Okamoto 2019)
			nx = obj.nx;
			nu = obj.nu;
			nw = obj.nw;
			N = obj.N;
			A_sys = obj.A;
			B_sys = obj.B;
			G_sys = obj.G; % G is the process noise input matrix
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
		
	end
end

