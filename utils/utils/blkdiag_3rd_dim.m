function out = blkdiag_3rd_dim(x)
    % blkdiag_3rd_dim(x) - unpacks 3rd dimension of x and calls blkdiag on
    % the 2d matrices. Can be used for cvx variables.
    c = num2cell(x, [1 2]);
    out = blkdiag(c{:}); 
end