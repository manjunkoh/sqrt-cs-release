function [K, P] = solve_lqr_tv(A, B, Q, R, Q_f, options)
%SOLVE_LQR_TV Solve finite-horizon time-varying LQR via Riccati equations
%
%   [K, P] = solve_lqr_tv(A, B, Q, R, Q_f)
%   [K, P] = solve_lqr_tv(A, B, Q, R, Q_f, options)
%
%   Solves the discrete-time finite-horizon time-varying LQR problem:
%   
%   minimize J = x_N^T Q_f x_N + sum_{k=0}^{N-1} [x_k^T Q_k x_k + u_k^T R_k u_k]
%   subject to x_{k+1} = A_k x_k + B_k u_k
%
%   The optimal control law is: u_k = -K_k x_k
%
%   Inputs:
%       A       - System dynamics matrix (nx x nx x N) or (nx x nx) for time-invariant
%       B       - Input matrix (nx x nu x N) or (nx x nu) for time-invariant
%       Q       - State cost matrix (nx x nx x N) or (nx x nx) for time-invariant
%       R       - Control cost matrix (nu x nu x N) or (nu x nu) for time-invariant
%       Q_f     - Terminal state cost matrix (nx x nx) [optional, default: Q]
%       options - Struct with optional fields:
%           .verbose - Display progress (default: false)
%
%   Outputs:
%       K       - Feedback gain matrices (nu x nx x N)
%       P       - Riccati matrices (nx x nx x N+1), with P(:,:,N+1) = Q_f
%
%   Example:
%       % Time-invariant system
%       A = [1 1; 0 1]; B = [0.5; 1];
%       Q = eye(2); R = 1; Q_f = 10*eye(2);
%       [K, P] = solve_lqr_tv(A, B, Q, R, Q_f);
%
%       % Time-varying system
%       nx = 2; nu = 1; N = 10;
%       A = repmat([1 1; 0 1], [1, 1, N]);
%       B = repmat([0.5; 1], [1, 1, N]);
%       Q = repmat(eye(2), [1, 1, N]);
%       R = repmat(1, [1, 1, N]);
%       [K, P] = solve_lqr_tv(A, B, Q, R, Q_f);

    arguments
        A (:,:,:) double
        B (:,:,:) double
        Q (:,:,:) double
        R (:,:,:) double
        Q_f (:,:) double = []
        options.verbose (1,1) logical = false
    end
    
    % Get dimensions
    [nx, ~, N_A] = size(A);
    [~, nu, N_B] = size(B);
    [~, ~, N_Q] = size(Q);
    [~, ~, N_R] = size(R);
    
    % Determine time horizon N
    N = max([N_A, N_B, N_Q, N_R]);
    
    % Check if matrices are time-invariant (size 2) and expand to 3D if needed
    if size(A, 3) == 1
        A = repmat(A, [1, 1, N]);
    end
    if size(B, 3) == 1
        B = repmat(B, [1, 1, N]);
    end
    if size(Q, 3) == 1
        Q = repmat(Q, [1, 1, N]);
    end
    if size(R, 3) == 1
        R = repmat(R, [1, 1, N]);
    end
    
    % Set default terminal cost if not provided
    if isempty(Q_f)
        Q_f = Q(:,:,end);
    end
    
    % Validate dimensions
    if size(Q_f, 1) ~= nx || size(Q_f, 2) ~= nx
        error('solve_lqr_tv: Q_f must be %d x %d', nx, nx);
    end
    
    % Initialize output arrays
    P = zeros(nx, nx, N+1);
    K = zeros(nu, nx, N);
    
    % Terminal condition: P_N = Q_f
    P(:,:,N+1) = Q_f;
    
    % Backward Riccati recursion: k = N, N-1, ..., 1
    for k = N:-1:1
        A_k = A(:,:,k);
        B_k = B(:,:,k);
        Q_k = Q(:,:,k);
        R_k = R(:,:,k);
        P_kp1 = P(:,:,k+1);
        
        % Compute intermediate matrices for Riccati equation
        % P_k = Q_k + A_k^T P_{k+1} A_k - A_k^T P_{k+1} B_k (R_k + B_k^T P_{k+1} B_k)^{-1} B_k^T P_{k+1} A_k
        
        % Compute R + B^T P_{k+1} B
        R_BPB = R_k + B_k' * P_kp1 * B_k;
        
        % Check if R_BPB is invertible (should be positive definite)
        if min(real(eig(R_BPB))) < 1e-10
            warning('solve_lqr_tv: R + B^T P_{k+1} B near singular at k=%d. Adding regularization.', k);
            R_BPB = R_BPB + 1e-10 * eye(size(R_BPB, 1));
        end
        
        % Compute feedback gain: K_k = (R_k + B_k^T P_{k+1} B_k)^{-1} B_k^T P_{k+1} A_k
        K(:,:,k) = R_BPB \ (B_k' * P_kp1 * A_k);
        
        % Compute Riccati matrix: P_k = Q_k + A_k^T P_{k+1} (A_k - B_k K_k)
        P(:,:,k) = Q_k + (A_k - B_k * K(:,:,k))' * P_kp1 * (A_k - B_k * K(:,:,k)) + K(:,:,k)' * R_k * K(:,:,k);
        
        % Alternative (equivalent) form:
        % P(:,:,k) = Q_k + A_k' * P_kp1 * A_k - A_k' * P_kp1 * B_k * K(:,:,k);
        
        if options.verbose && mod(k, max(1, floor(N/10))) == 0
            fprintf('solve_lqr_tv: Completed k = %d / %d\n', k, N);
        end
    end
    
    if options.verbose
        fprintf('solve_lqr_tv: Completed backward recursion for N = %d time steps\n', N);
    end
end

