classdef ExampleClass4 < SCPProblem
    
    properties
        init_guess_struct
        impose_trust_region_struct = struct('x', true, 'u', false);
    end

	properties (SetAccess=immutable)
		x_init
		x_fin
	end

    properties (Constant)
        sigma = 10;  % Lorenz parameter
        rho = 28;    % Lorenz parameter
        beta = 8/3;  % Lorenz parameter
        t_N = 2;     % final time
        u_min = -50;  % minimum control input
        u_max = 50;   % maximum control input
        Nseg = 40;    % number of segments
    end

    properties
        DS
        DeltaT
        t_his
        A
        B
        c
    end

    methods

        function obj = ExampleClass4()
            % Must inherit from SCPProblem class
            obj@SCPProblem();

            obj.DS     = Lorenz(obj.sigma, obj.rho, obj.beta);
            obj.DeltaT = obj.t_N/(obj.Nseg);
            obj.t_his  = linspace(0, obj.t_N, obj.Nseg+1);

			% Set the initial and final states to the equilibrium points

			% Equilibrium points: origin (0,0,0) and symmetric points at (~±8.49, ±8.49, 27)
			% For standard parameters (σ=10, ρ=28, β=8/3):
			% - Equilibrium 1: (0, 0, 0) - unstable
			% - Equilibrium 2: (√(β(ρ-1)), √(β(ρ-1)), ρ-1) ≈ (8.49, 8.49, 27)
			% - Equilibrium 3: (-√(β(ρ-1)), -√(β(ρ-1)), ρ-1) ≈ (-8.49, -8.49, 27)
			obj.x_init = [sqrt(obj.beta * (obj.rho - 1)); 
						sqrt(obj.beta * (obj.rho - 1)); 
						obj.rho - 1];

			obj.x_fin = [-sqrt(obj.beta * (obj.rho - 1)); 
						-sqrt(obj.beta * (obj.rho - 1)); 
						obj.rho - 1];

			% Set the initial guess (linear interpolation between initial and final states)
			obj.init_guess_struct.x     = linspace_vec(obj.x_init, obj.x_fin, obj.Nseg+1);
			obj.init_guess_struct.u     = zeros(obj.DS.nu, obj.Nseg); % zero control initial guess

            % Must initialize the parent class
            obj.initialize()
        end

        function vars = define_vars(obj)
            vars.x = sdpvar(obj.DS.nx, obj.Nseg+1);
            vars.u = sdpvar(obj.DS.nu, obj.Nseg);
        end

        function J0 = objective(obj, vars)
            % Minimize control effort (integral of |u|)
            J0 = sum(abs(vars.u)) * obj.DeltaT;
        end

        function constraints = convex_eq(obj, vars)
            constraints = [
                vars.x(:,1)   == obj.x_init
                vars.x(:,end) == obj.x_fin
            ];
        end

        function constraints = convex_ineq(obj, vars)
            constraints = [
                obj.u_min <= vars.u
                vars.u <= obj.u_max
            ];
        end

        function update_parameters(obj, ref_vars)
            [obj.A, obj.B, obj.c] = obj.DS.discretize_LT(ref_vars.x, ref_vars.u, obj.t_his);
        end

        function constraintLHS = noncvx_eq(obj, vars)
            constraintLHS = [];
            % Dynamics equation
            for k = 1:obj.Nseg
                [~, x] = obj.DS.propagate_with_LT(vars.x(:,k), obj.t_his(k:k+1), vars.u(:,k));
                x_kplus1 = x(end,:)';
                constraintLHS = [constraintLHS; vars.x(:,k+1) - x_kplus1];
            end
        end

        function constraintLHS = noncvx_ineq(obj, vars)
            constraintLHS = [];
            % No non-convex inequality constraints for this problem
        end

        function constraintLHS = noncvx_eq_relaxed(obj, vars, ref_vars)
            constraintLHS = [];

            for k = 1:obj.Nseg
                constraintLHS = [constraintLHS
                vars.x(:,k+1) - obj.A(:,:,k) * vars.x(:,k) - obj.B(:,:,k)* vars.u(:,k) - obj.c(:,:,k)];
            end
        end

        function constraintLHS = noncvx_ineq_relaxed(obj, vars, ref_vars)
            constraintLHS = [];
            % No non-convex inequality constraints for this problem
        end

    end
end

