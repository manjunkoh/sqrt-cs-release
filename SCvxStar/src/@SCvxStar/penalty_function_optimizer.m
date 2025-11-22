% Called when creating the optimizer objective
function P = penalty_function_optimizer(obj, xi, zeta, lambda, mu, w, xi_norm_squared, zeta_norm_squared)
    switch obj.constParams.penalty_method
        case {'AL'}
            P = dot(lambda, xi) + dot(mu, zeta)+  w/2 * xi_norm_squared + w/2 * zeta_norm_squared;
    otherwise
        error('penalty_method except AL not supported for optimizer');
    end
end