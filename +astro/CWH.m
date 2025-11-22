%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CWH (Clohessy-Wiltshire-Hill) equation class definition
% Requires the DynamicalSystem class 
% 
% The CWH equations model relative orbital motion in proximity to a 
% chief satellite in a circular orbit. The state is x = [r; v] where
% r and v are 3D position and velocity vectors in the rotating frame.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef CWH < astro.DynamicalSystem

    properties (SetAccess=immutable)
        n % Mean motion: n = sqrt(μ/r_0^3)
        mu % Gravitational parameter (default: Earth)
        A % A matrix
        B % B matrix
    end
    properties (Constant)
        is_control_affine = true;
        nx = 6;  % State: [r; v] where r and v are 3D
        nu = 3;  % Control: acceleration in 3D
        B_impulse = [];
    end

    methods
        function obj = CWH(r0, mu)
            % CWH Constructor
            % Inputs:
            %   r0: Chief orbit radius (km or m, must match mu units)
            %   mu: Gravitational parameter (optional, defaults to Earth in km^3/s^2)
            arguments
                r0
                mu = 3.986004418e5  % Earth's gravitational parameter in km^3/s^2
            end
            obj.mu = mu;
            obj.n = sqrt(mu / r0^3);
            obj.A = obj.get_A_matrix();
            obj.B = obj.get_B_matrix();
        end   
        
        function A = get_A_matrix(obj)
            % Get the A matrix for CWH dynamics
            % A = [[0_3x3, I_3], [A_21, A_22]]
            n = obj.n;
            A_21 = [3*n^2, 0, 0;
                    0, 0, 0;
                    0, 0, -n^2];
            A_22 = [0, 2*n, 0;
                    -2*n, 0, 0;
                    0, 0, 0];
            A = [zeros(3,3), eye(3);
                 A_21, A_22];
        end
        
        function B = get_B_matrix(obj)
            % Get the B matrix for CWH dynamics (control input matrix)
            % Control is acceleration, so B = [0_3x3; I_3]
            B = [zeros(3,3); eye(3)];
        end

        function dxdt = EOM(obj, t, x)
            % Equations of motion: dx/dt = Ax + Bu
            % For EOM without control input, u = 0
            dxdt = obj.A * x;
        end

        function dfdx = dfdx(obj, t, x)
            % Jacobian with respect to state (constant for CWH)
            dfdx = obj.A;
        end

        function dfdu = dfdu(obj, t, x)
            % Jacobian with respect to control (constant for CWH)
            dfdu = obj.B;
        end

        function add_impulse(obj)
            error('Impulsive control not implemented for CWH')
        end

    end
end

