function out = quantile_from_gx2(x, mu, P)
    % Compute quantile of the norm of a Gaussian random vector
    % x: quantile
    % mu: mean
    % P: covariance
    % out = gx2cdf(p, reshape(sqrt(diag(P)), 1, []), 1, norm(x)^2, 0, 0);

    % gx2cdf(x,w,k,lambda,m,s)
    % x         points at which to evaluate the CDF
	% w         row vector of weights of the non-central chi-squares
	% k         row vector of degrees of freedom of the non-central chi-squares
	% lambda    row vector of non-centrality paramaters (sum of squares of
	%           means) of the non-central chi-squares
	% m         mean of normal term
    % s         sd of normal term
    
    k = 1;
    m = 0;
    s = 0;

    % w = reshape(diag(P), 1, []);
    % lambda = reshape(mu.^2./diag(P), 1, []);
 
    [V, D] = eig(P);
    w = diag(D)';
    lambda = ((V'*inv(sqrtm(P))*mu).^2)';

    out = sqrt(gx2cdf(x,w,k,lambda,m,s));
