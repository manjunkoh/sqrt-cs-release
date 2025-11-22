%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DynamicalSystem is a class that defines the interface for dynamical systems.
% It is an abstract class that must be subclassed to be used.
% Child classes must implement the EOM, dfdx, dfdu, and add_impulse methods.
% The class provides methods for propagating the system, computing the
% state transition matrix, and discretizing the system for optimization.
% The current implementation assumes that the system is control affine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef DynamicalSystem < dynamicprops 
    
    properties
        odeopts struct = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
        integrator function_handle = @ode78;
        SF astro.ScaleFactor
    end

    properties (Abstract, Constant)
        nx
        nu
        B_impulse
    end

    methods 

        function obj = DynamicalSystem(SF)
            if nargin < 1
                SF = astro.ScaleFactor(1,1);
            end
            obj.SF = SF;        
        end

        function varargout = propagate(obj, x0, tspan, odeopts)
            arguments
                obj
                x0 
                tspan
                odeopts = obj.odeopts
            end

            [varargout{1:nargout}] = obj.integrator(@obj.EOM, tspan, x0, odeopts);
        end

        function y0 = get_initial_state(obj, x0)
            switch size(x0, 1)
                case obj.nx
                    y0 = [x0; reshape(eye(obj.nx), [], 1)];
                case obj.nx^2 + obj.nx
                    y0 = x0;
                otherwise
                    error('Initial condition must be either state or augmented state.');
            end
        end

        function varargout = propagate_with_STM(obj, x0, tspan, odeopts)
            arguments
                obj
                x0 
                tspan
                odeopts = obj.odeopts
            end
            y0 = obj.get_initial_state(x0);
            [varargout{1:nargout}] = obj.integrator(@obj.EOM_with_STM, tspan, y0, odeopts);
        end

        function varargout = propagate_with_LT(obj, x0, tspan, u, odeopts)
            arguments
                obj
                x0
                tspan
                u
                odeopts = obj.odeopts
            end

            if isa(u, "double")
                [varargout{1:nargout}] = obj.integrator(@(t, x) obj.EOM_with_LT(t, x, u), tspan, x0, odeopts);
            elseif isa(u, "function_handle")
                [varargout{1:nargout}] = obj.integrator(@(t, x) obj.EOM_with_LT(t, x, u(t)), tspan, x0, odeopts);
            end
        end

        function varargout = propagate_with_LT_with_STM(obj, x0, tspan, u, odeopts)
            arguments
                obj
                x0
                tspan
                u
                odeopts = obj.odeopts
            end

            y0 = obj.get_initial_state(x0);
            % TODO: include this if statement for all functions
            if isa(u, "double")
                [varargout{1:nargout}] = obj.integrator(@(t, y) obj.EOM_with_LT_with_STM(t, y, u), tspan, y0, odeopts);
            elseif isa(u, "function_handle")
                [varargout{1:nargout}] = obj.integrator(@(t, y) obj.EOM_with_LT_with_STM(t, y, u(t)), tspan, y0, odeopts);
            end
        end

        function varargout = propagate_with_STM_with_covariance(obj, x0, P0, tspan, Q, odeopts)
            arguments
                obj
                x0
                P0
                tspan
                Q
                odeopts = obj.odeopts
            end

            y0 = [x0; reshape(eye(obj.nx), [], 1); reshape(P0, [], 1)];
            [varargout{1:nargout}] = obj.integrator(@(t, y) obj.EOM_with_STM_with_covariance(t, y, Q), tspan, y0, odeopts);
        end

        function [dydt, A] = EOM_with_STM(obj, t, y)
            x = y(1:obj.nx);
            STM = reshape(y(obj.nx+1:end), obj.nx, obj.nx);
            A = obj.dfdx(t, x);
            dydt = [obj.EOM(t, x); reshape(A*STM, [], 1)];
        end

        function dxdt = EOM_with_LT(obj, t, x, u)
            dxdt = obj.EOM(t, x) + obj.dfdu(t, x)*u;
        end

        function dydt = EOM_with_LT_with_STM(obj, t, y, u)
            x = y(1:obj.nx);
            STM = reshape(y(obj.nx+1:end), obj.nx, obj.nx);

            A = obj.dfdx(t, x);
            B = obj.dfdu(t, x);
            dydt = [obj.EOM(t, x) + B*u; reshape(A*STM, [], 1)];
        end

        function dydt = EOM_with_STM_with_covariance(obj, t,y,Q)
            x = y(1:obj.nx);
            STM = reshape(y(obj.nx+1:obj.nx+obj.nx^2), obj.nx, obj.nx);
            P = reshape(y(obj.nx+1+obj.nx^2:end), obj.nx, obj.nx);

            A = obj.dfdx(t, x);
            dydt = [obj.EOM(t, x); reshape(A*STM, [], 1); reshape(A*P + P*A' + Q, [], 1)];
        end
            
        % Methods in separate files
        % If there is an error 'Method is not defined or is removed from MATLAB search path',
        % this is a known bug: 
        % https://www.mathworks.com/matlabcentral/answers/467384-method-is-not-defined-or-is-removed-from-matlab-search-path
        % Restarting MATLAB should fix the issue.

        [A, B, c, G, trajs, GG] = discretize_LT_SDE(obj, x_ref, u_ref, t_his, g)
        [A_, B_, c_, G_, trajs] = discretize_impulsive_SDE(obj, x_ref, u_ref, t_his, g)

    end

    methods (Abstract)
        dxdt = EOM(obj, t, x)
        A = dfdx(obj, t, x)
        B = dfdu(obj, t, x)
        x_plus = add_impulse(obj, x, u)
    end

end