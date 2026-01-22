function fill_confidence_ellipse(mu, P, percentile, varargin)
    
	[V,D] = eig(P(1:2,1:2));

	s = sqrt(chi2inv(percentile, 2));
	a = s * sqrt(D(1,1));
	b = s * sqrt(D(2,2));

	theta = atan2(V(2,1), V(1,1));

	fillEllipse(a, b, mu(1:2), theta, varargin{:});
end

function fillEllipse(a,b,c,theta,varargin)
    t = linspace(0, 2 * pi);

    x = c(1) + a*cos(t)*cos(theta) - b*sin(t)*sin(theta);
    y = c(2) + a*cos(t)*sin(theta) + b*sin(t)*cos(theta);

    fill(x,y, varargin{:});
end