%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lorenz system class definition
% Requires the DynamicalSystem class 
% The Lorenz system with control input u added to the dot{y} term:
%   dx/dt = sigma * (y - x)
%   dy/dt = rho * x - y - x*z + u
%   dz/dt = x*y - beta * z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Lorenz < astro.DynamicalSystem

    properties (SetAccess=immutable)
        sigma % parameter sigma
        rho   % parameter rho
        beta  % parameter beta
    end
    properties (Constant)
        is_control_affine = true;
        nx = 3;
        nu = 1;
        B_impulse = [];
    end

    methods
        function obj = Lorenz(sigma, rho, beta)
            if nargin < 1
                sigma = 10; % default value
            end
            if nargin < 2
                rho = 28; % default value
            end
            if nargin < 3
                beta = 8/3; % default value
            end
            obj.sigma = sigma;
            obj.rho = rho;
            obj.beta = beta;
        end   

        function dxdt = EOM(obj, t, x)
            % Equations of motion without control
            x_val = x(1);
            y_val = x(2);
            z_val = x(3);
            dxdt = [obj.sigma * (y_val - x_val);
                    obj.rho * x_val - y_val - x_val * z_val;
                    x_val * y_val - obj.beta * z_val];
        end

        function dfdx = dfdx(obj, t, x)
            % Jacobian with respect to state
            x_val = x(1);
            y_val = x(2);
            z_val = x(3);
            dfdx = [-obj.sigma,          obj.sigma,           0;
                    obj.rho - z_val,     -1,                  -x_val;
                    y_val,                x_val,               -obj.beta];
        end

        function dfdu = dfdu(obj, t, x)
            % Jacobian with respect to control (input to dot{y} term)
            dfdu = [0; 1; 0];
        end

        function x_plus = add_impulse(obj, x, u)
            error('Impulsive control not implemented for Lorenz system')
        end

    end
end

