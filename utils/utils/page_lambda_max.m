function lmd_arr = page_lambda_max(X)
    lmd_arr = NaN(size(X,3),1);
    for i = 1:size(X,3)
        lmd_arr(i) = lambda_max(X(:,:,i));
    end
end