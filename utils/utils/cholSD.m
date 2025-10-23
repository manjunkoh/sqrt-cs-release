function L = cholSD(M)
    try
        L = chol(M, 'lower');
    catch
        [V, D] = eig(M);
        L = V * sqrt(D);
    end
end