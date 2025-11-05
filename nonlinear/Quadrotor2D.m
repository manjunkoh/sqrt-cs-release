%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quadrotor2D class definition
% Requires the DynamicalSystem class 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Quadrotor2D < astro.DynamicalSystem

    properties (SetAccess=immutable)
        k_D % drag coefficient
    end
    properties (Constant)
        is_control_affine = true;
        nx = 4;
        nu = 2;
        B_impulse = [];
    end

    methods
        function obj = Quadrotor2D(k_D)
            obj.k_D = k_D;
        end   

        function dxdt = EOM(obj, t, x)
            v = x(3:4);
            dxdt = [v; - obj.k_D * norm(v) * v];
        end

        function dfdx = dfdx(obj, t, x)
            v = x(3:4);
            v_norm = norm(v);
            dfdx = [zeros(2,2),  eye(2);
                    zeros(2,2),  (-obj.k_D * (v_norm * eye(2) + (v * v') / v_norm))];
        end

        function dfdu = dfdu(obj, t, x)
            dfdu = [zeros(2,2); eye(2)];
        end

        function add_impulse(obj)
            error('Impulsive control not implemented for Quadrotor2D')
        end

    end
end