classdef CovarianceSteeringBase < handle
	% Base class for Covariance Steering optimization problems
	% 
	% This is an abstract base class that defines common properties and methods
	% for different covariance steering solution methods.
	%
	% Derived classes must implement:
	%   - set_sdpvars()
	%   - set_objective()
	%   - set_constraints()
	%   - set_feedback_gains()
	
	properties
		A % System dynamics matrix
		B % Input matrix
		G % Process noise input matrix (Cholesky factor of process noise covariance)
		P_0 % Initial state covariance
		P_f % Final state covariance
		P_ref % Reference covariance for linearization (e.g., from previous SCP iteration), size (nx x nx x N+1) or empty; only used for FullCovarianceSteering
		Q % State cost matrix; assume constant over time
		R % Control cost matrix; assume constant over time
		nx % State dimension
		nu % Control dimension
		nw % Process noise dimension
		N % Time horizon
		sdpvars % Struct to hold YALMIP sdpvar variables
		covariance_scaling = 1; % Scaling factor for state covariance
		K = []; % Feedback gains; computed after solving
		optimal_objective = NaN; % Optimal objective value after solving
		% Chance constraint support (optional, for derived classes that implement it)
		% State constraints: cell array of structs with fields:
		%   - type: 'affine' or 'norm'
		%   - For affine: alpha (vector), beta (scalar), p (violation prob), nodes (optional, default: all)
		%   - For norm: gamma (scalar), p (violation prob), n (dimension), nodes (optional, default: all)
		chance_constraints_state = {}
		mu_0 = []  % initial mean (nx x 1) or empty -> assumed zero
		mu_f = []  % terminal mean (nx x 1) or empty -> assumed zero
		waypoints = {}  % cell array of waypoint structs with fields: 'node' (scalar, 1:N+1) and 'mu' (nx x 1 vector)
		mu % state mean trajectory, set after solving (if applicable)
		v  % control mean trajectory, set after solving (if applicable)
		P  % state covariance matrices, set after solving (nx x nx x N+1)
		P_u % control covariance matrices, set after solving (nu x nu x N)
	end
	
	methods
		function obj = CovarianceSteeringBase(options)
			arguments
				options.A
				options.B
				options.G
				options.P_0
				options.P_f
				options.Q
				options.R
				options.N
				options.covariance_scaling = 1
				options.chance_constraints_state = {}
				options.mu_0 = []
				options.mu_f = []
				options.waypoints = {}
				options.P_ref = []
			end
			
			obj.A = options.A;
			obj.B = options.B;
			obj.G = options.G;
			obj.P_0 = options.P_0;
			obj.P_f = options.P_f;
			obj.P_ref = options.P_ref;
			obj.Q = options.Q;
			obj.R = options.R;
			obj.nx = size(options.A, 1);
			obj.nu = size(options.B, 2);
			obj.nw = size(options.G, 2);
			obj.N = options.N;
			obj.covariance_scaling = options.covariance_scaling;
			
			% Optional chance constraint and mean parameters
			obj.chance_constraints_state = options.chance_constraints_state;
			if ~isempty(options.mu_0)
				obj.mu_0 = options.mu_0(:); % ensure column vector
			end
			if ~isempty(options.mu_f)
				obj.mu_f = options.mu_f(:); % ensure column vector
			end
			obj.waypoints = options.waypoints;
		end

		function diagnostic = solve(obj, options)
			% Solve the SDP using YALMIP
			% This method calls the abstract methods implemented by derived classes
			arguments
				obj
				options.verbose = 0
				options.solver = 'mosek'
				options.savesolveroutput = true
			end
			
			% Call abstract methods implemented by derived classes
			obj.set_sdpvars();
			obj.set_objective();
			obj.set_constraints();

			settings = sdpsettings();
			settings.verbose = options.verbose;
			settings.solver = options.solver;
			settings.savesolveroutput = options.savesolveroutput;

			diagnostic = optimize(obj.sdpvars.constraints, obj.sdpvars.J, settings);

			switch diagnostic.problem
				case 0
					obj.set_feedback_gains();
					obj.optimal_objective = value(obj.sdpvars.J);
				case 1
					warning('YALMIP:Infeasible', 'Problem infeasible');
				case 4
					obj.set_feedback_gains();
					obj.optimal_objective = value(obj.sdpvars.J);
				otherwise
					warning('YALMIP:Failed', 'Solver failed with code %d', diagnostic.problem);
			end
		end

		function varargout = plot_Y_lambda_max(obj, varargin)
			N = obj.N;
			
			if isempty(obj.P_u)
				error('CovarianceSteeringBase: P_u not set. Solve the problem first.');
			end
			
			% Compute maximum eigenvalues of control covariance matrices
			Y_lmd_max = NaN(N, 1);
			for k = 1:N
				% Ensure matrix is symmetric for eigenvalue computation
				P_u_k = (obj.P_u(:,:,k) + obj.P_u(:,:,k)') / 2;
				Y_lmd_max(k) = lambda_max(P_u_k);
			end
			
			[varargout{1:nargout}] = stairsZOH(0:N, Y_lmd_max, varargin{:});
			xlabel('Node $k$')
			ylabel('$\lambda_{\max}(Y_k)$')
		end

		function plot_state_sigmas(obj, varargin)
			nx = obj.nx;
			N = obj.N;
			
			if isempty(obj.P)
				error('CovarianceSteeringBase: P not set. Solve the problem first.');
			end
			
			% Compute standard deviations from covariance matrices
			sigmas = NaN(nx, N+1);
			for k = 1:N+1
				for i = 1:nx
					sigmas(i,k) = sqrt(obj.P(i,i,k));
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

		function J_step = get_objective(obj)
			% Get objective function value per step
			nx = obj.nx;
			nu = obj.nu;
			N = obj.N;
			
			if isempty(obj.P) || isempty(obj.P_u)
				error('CovarianceSteeringBase: P or P_u not set. Solve the problem first.');
			end
			
			J_step = zeros(N, 1);
			
			for k = 1:N
				state_cost = trace(obj.Q * obj.P(:,:,k));
				control_cost = trace(obj.R * obj.P_u(:,:,k));
				J_step(k) = state_cost + control_cost;
			end
		end
		
		% Abstract methods that must be implemented by derived classes
		set_sdpvars(obj)
		set_objective(obj)
		set_constraints(obj)
		set_feedback_gains(obj)
	end
end