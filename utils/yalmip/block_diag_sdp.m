function L = block_diag_sdp(x,y,N)
% returns a (Nx,Ny) sdpvar of a (x,y) block N by N diagonal matrix
% more efficient than naively indexing the matrix to set 0 constraints
L_block = sdpvar(repmat(x,1,N), repmat(y,1,N));
L = blkdiag(L_block{:});
end