function arr = lambda_maxes(X)
    arr = NaN(size(X,3),1);
    for i = 1:size(X,3)
        arr(i) = lambda_max(X(:,:,i));
    end
end