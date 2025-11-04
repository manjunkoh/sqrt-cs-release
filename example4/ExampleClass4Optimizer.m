classdef ExampleClass4Optimizer < ExampleClass4
    
    properties 
        changing_params = struct();
        is_used_in_constraint_update = struct('x', true, 'u', false);
    end

    methods

        function obj = ExampleClass4Optimizer()
            % Must inherit from ExampleClass4
            obj@ExampleClass4();
        end

        function val = get_use_optimizer(obj)
            val = true;
        end

        function set_changing_parameters_sdp(obj)
            obj.changing_params.A = sdpvar(obj.DS.nx, obj.DS.nx, obj.Nseg, 'full');
            obj.changing_params.B = sdpvar(obj.DS.nx, obj.DS.nu, obj.Nseg, 'full');
            obj.changing_params.c = sdpvar(obj.DS.nx, 1, obj.Nseg);
        end

        function p = get_changing_parameters(obj, ref_vars)
            % Numerical update of the parameters that appear in the constraints
            [p.A, p.B, p.c] = obj.DS.discretize_LT(ref_vars.x, ref_vars.u, obj.t_his);
        end

        function constraints = get_changing_constraints(obj, vars, ref_vars)
            obj.set_changing_parameters_sdp();
            p = obj.changing_params;

            % Vectorized version of the dynamics constraints
            nx = obj.DS.nx;
            nu = obj.DS.nu;
            Nseg = obj.Nseg;

            A_blk = [eye(nx), zeros(nx, nx*Nseg); blkdiag_along_3rd_dim(p.A), zeros(nx*Nseg, nx)];
            B_blk = [zeros(nx, nu*Nseg); blkdiag_along_3rd_dim(p.B)];
            c_blk = [zeros(nx, 1); p.c(:)];
            
            constraints_eom = (eye(nx*(Nseg+1)) - A_blk) * vars.x(:) - B_blk * vars.u(:) - c_blk;
            constraints_eom = constraints_eom(nx+1:end);

            constraints = [constraints_eom == obj.slack_noncvx_eq]; 
        end

    end
end

function out = blkdiag_along_3rd_dim(x)
    % unpacks 3rd dimension of x and calls blkdiag on the 2d matrices. 
    c = num2cell(x, [1 2]);
    out = blkdiag(c{:}); 
end

