function L = page_tril_sdp(n,N)
% returns a (n,n,N) sdpvar of a where each page is a (n,n) lower triangular matrix
L = repmat(tril(ones(n)), [1 1 N]);
L = L .* sdpvar(n,n,N,'full');
end