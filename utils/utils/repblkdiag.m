function repM = repblkdiag(M, n)
    C = repmat({M},n,1);
    repM = blkdiag(C{:});
end
