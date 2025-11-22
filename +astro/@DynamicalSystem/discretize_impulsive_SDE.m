% Note: the state x_k is the state BEFORE the impulse u_k is applied
function [A_, B_, c_, G_, trajs] = discretize_impulsive_SDE(obj, x_ref, u_ref, t_his, g)

    arguments
        obj
        x_ref
        u_ref
        t_his
        g function_handle 
    end

    nx = obj.nx;
    nu = obj.nu;
    N = length(t_his) - 1;

    A_ = NaN(nx,nx,N);
    B_ = NaN(nx,nu,N);
    c_ = NaN(nx,1,N);
    G_ = NaN(nx,nx,N);
    trajs = cell(N,1);

    for k = 1:N
        [A_(:,:,k), c_(:,:,k), G_(:,:,k), trajs{k}] = discretize_segment(obj, x_ref(:,k), u_ref(:,k), t_his(k:k+1), g);
        B_(:,:,k) = A_(:,:,k) * obj.B_impulse;
    end
    % B_ = repmat(obj.B_impulse, [1 1 N+1]);

end

function [A, c, G, traj] = discretize_segment(obj, x_ref_k, u_ref_k, tspan, g)
    nx = obj.nx;
    
    Y = [x_ref_k + obj.B_impulse * u_ref_k
        vec(eye(nx))
        zeros(nx, 1)
        zeros(nx^2, 1)
        ];

    traj = obj.integrator(@(t, y) EOM_with_discrete_matrices(obj, t, y, g), ...
                            tspan, ...
                            Y, ...
                            obj.odeopts);
    y_f = traj.y(:,end);
    [A, c, G] = unpack_matrices(y_f, nx);
end

function dYdt = EOM_with_discrete_matrices(obj, t, Y, g)

    x = Y(1:obj.nx);
    STM = reshape(Y(obj.nx + 1:obj.nx + obj.nx ^ 2), [obj.nx, obj.nx]);

    A = obj.dfdx(t, x);
    c = obj.EOM(t, x) - A * x;
    G = g(t, x);

    dYdt = [
            obj.EOM(t, x);
            vec(A * STM)
            STM \ c
            vec((STM\G) * (STM\G)')
            ];
end
    
function [A, c, G] = unpack_matrices(y, nx)
    A = reshape(y(nx+1:nx+nx^2), [nx nx]);
    c = A * y(nx+nx^2+1:2*nx+nx^2);
    G = A * chol(reshape(y(2*nx+nx^2+1:end), [nx nx]), 'lower');
end