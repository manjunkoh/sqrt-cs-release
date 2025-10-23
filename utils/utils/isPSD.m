function out = isPSD(M)
    out = 1;
    try 
        chol(M);
    catch
        out = 0;
    end
end