function C = cholAll(M, varargin)

    C = zeros(size(M));
    for i = 1:size(M,3)
        try
            C(:,:,i) = chol(M(:,:,i), varargin{:});
        catch
            C(:,:,i) = NaN(size(M(:,:,i)));
        end        
    end
end