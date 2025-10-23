function [out, idx_notPSD] = isPSD_allpages(X)
    D = pageeig(X);
    out = all(D(:) >= 0);
    
    if nargout > 1
        idx_notPSD = [];
        for i = 1:size(X,3)
            if any(D(:,i) < 0)
                idx_notPSD = [idx_notPSD i];
            end
        end
    end
end
