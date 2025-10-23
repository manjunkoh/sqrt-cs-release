function out = sum_trace(X)
    out = 0;
    for i = 1:size(X,3)
        out = out + trace(X(:,:,i));
    end
end