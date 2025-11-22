function [A, B, c, G, trajs, GG] = discretize_LT_SDE(obj, x_ref, u_ref, t_his, g)

    nx = obj.nx;
    nu = obj.nu;

    N = length(t_his) - 1;

    A = NaN(nx,nx,N);
    B = NaN(nx, nu, N);
    c = NaN(nx,1,N);
    G = NaN(nx, nx, N);
    GG = NaN(nx, nx, N);
    trajs = cell(N,1);

    for k = 1:N
        [A(:,:,k), B(:,:,k), c(:,:,k), G(:,:,k), trajs{k}, GG(:,:,k)] = discretize_segment(obj, x_ref(:,k), u_ref(:,k), t_his(k:k+1), g);
    end

end

function [A_k, B_k, c_k, G_k, traj, GG_k] = discretize_segment(obj, x_ref_k, u_ref_k, tspan, g)
    nx = obj.nx;
    nu = obj.nu;
    
    Y_k = [x_ref_k
        vec(eye(nx))
        zeros(nx * nu, 1)
        zeros(nx^2, 1)
        ];

    traj = obj.integrator(@(t, y) EOM_with_LT_with_discrete_matrices(obj, t, y, u_ref_k, g), ...
                            tspan, ...
                            Y_k, ...
                            obj.odeopts);
    y_kp1 = traj.y(:,end);
    x_kp1 = y_kp1(1:nx);
    A_k = reshape(y_kp1(nx+1:nx+nx^2), [nx nx]);
    B_k = reshape(y_kp1(nx+nx^2+1:nx+nx^2+nx*nu), [nx nu]);
    c_k = x_kp1 - A_k * x_ref_k - B_k * u_ref_k;
    GG_k = reshape(y_kp1(nx+nx^2+nx*nu+1:end), [nx nx]);
    G_k = chol(GG_k, 'lower');

end

function dYdt = EOM_with_LT_with_discrete_matrices(obj, t, Y, u_k, g)

    x = Y(1:obj.nx);
    Phi_A = reshape(Y(obj.nx + 1:obj.nx + obj.nx ^ 2), [obj.nx, obj.nx]);
    Phi_B = reshape(Y(obj.nx + obj.nx ^ 2 + 1:obj.nx + obj.nx ^ 2 + obj.nx * obj.nu), [obj.nx, obj.nu]);
    Phi_Q = reshape(Y(obj.nx + obj.nx ^ 2 + obj.nx * obj.nu + 1:end), [obj.nx, obj.nx]);

    dfdx = obj.dfdx(t, x);
    dfdu = obj.dfdu(t, x);
    G = g(t, x);

    dYdt = [
            obj.EOM_with_LT(t, x, u_k)
            vec(dfdx * Phi_A)
            vec(dfdx * Phi_B + dfdu)
            vec(dfdx * Phi_Q + Phi_Q * dfdx.' + G * G.')
        ];
end