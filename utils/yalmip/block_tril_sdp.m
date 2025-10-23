function L = block_tril_sdp(x,y,N)
% returns a (Nx,Ny) sdpvar of a (x,y) block N by N lower triangular matrix
% more efficient than naively indexing the matrix to set 0 constraints
% https://groups.google.com/g/yalmip/c/PHnGnHXNcEA
L = kron(tril(ones(N)),ones(x,y));
L = L.*sdpvar(N*x,N*y);
end