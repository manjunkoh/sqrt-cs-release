function plot_confidence_ellipse(mu, P, percentile, varargin)
    
	[V,D] = eig(P(1:2,1:2));

	s = sqrt(chi2inv(percentile, 2));
	a = s * sqrt(D(1,1));
	b = s * sqrt(D(2,2));

	theta = atan2(V(2,1), V(1,1));

	plotEllipse(a, b, mu(1:2), theta, varargin{:});
end