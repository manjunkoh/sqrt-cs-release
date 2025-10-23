function [Qx, Rx, Qz, Rz, dQ, dR] = d_qr(X, dX)
%D_QR Computes the differentials of the QR decomposition.
%   This is a MATLAB implementation of the 'd_qr' R function from the paper
%   "Differentiating the QR Decomposition" by Jan de Leeuw.

    [n, m] = size(X);
    Z = X + dX;

    % Perform economy-size QR decomposition of the original matrix X and
    % the perturbed matrix Z = X + dX.
    [Qx_raw, Rx_raw] = qr(X, "econ");
    [Qz_raw, Rz_raw] = qr(Z, "econ");

    % MATLAB's economy qr(A) returns an R that is size(A). We need R to be
    % a square m x m matrix.
    % TODO: Is this step necessary? It seems redundant.
    Rx_temp = Rx_raw(1:m, :);
    Rz_temp = Rz_raw(1:m, :);

    % For uniqueness and consistent comparison, enforce the convention that R
    % has positive diagonal elements. This fixes sign ambiguities in Q and R.
    s_X = diag(sign(diag(Rx_temp)));
    Qx = Qx_raw * s_X;
    Rx = s_X * Rx_temp;

    s_Z = diag(sign(diag(Rz_temp)));
    Qz = Qz_raw * s_Z;
    Rz = s_Z * Rz_temp;

    % Find an orthonormal basis for the orthogonal complement of the column
    % space of Qx. This is a direct translation of the R code:
    % qp <- qr.Q(qr(cbind(qx, diag(n))))[, -(1:m)]
    [Q_full, ~] = qr([Qx, eye(n)]);
    Qp = Q_full(:, m+1:end);

    % Calculate the inverse of Rx. Using backslash is numerically stabler.
    Rxinv = Rx \ eye(m);

    % Calculate the intermediate matrix V
    V = Qx' * dX * Rxinv;

    % Create the anti-symmetric matrix 'A' from the lower triangular part of V.
    % The paper's lt(V) function is equivalent to MATLAB's tril(V).
    % A <- lt(V) - t(lt(V))
    M = tril(V);
    A = M - M';

    dR = (V - A) * Rx;

    % Calculate matrix 'B'
    B = Qp' * dX * Rxinv;
    dQ = Qx * A + Qp * B;
end