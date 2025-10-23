function x_interp = linspace_mat(X1, X2, n)
    assert(all(size(X1) == size(X2)), 'X1 and X2 must have the same dimensions');
    x_interp = NaN(size(X1, 1), size(X1, 2), n);
    for i = 1:size(X1, 1)
        for j = 1:size(X1, 2)
            x_interp(i,j,:) = linspace(X1(i,j), X2(i,j), n);
        end
    end
end