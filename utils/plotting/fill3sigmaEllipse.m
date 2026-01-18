function varargout = fill3sigmaEllipse(mu, cov, varargin)
%PLOT3SIGMAELLIPSE plots a 3 sigma ellipse

	assert(all(eig(cov(1:2, 1:2))));
	[V,d] = eig(cov(1:2, 1:2));

	V1 = V(:,1);
	try
		theta = atan2(V1(2), V1(1)); % angle of the largest eigenvector with the x-axis
	catch
		warning('Small imaginary value in covariance matrix eigenvalue.');
		return
	end

	s = 3;
	a = s * sqrt(d(1,1));
	b = s * sqrt(d(2,2));

	[varargout{1:nargout}] = fillEllipse(a,b,mu(1:2), theta, varargin{:});

end

function varargout = fillEllipse(a,b,c,theta,varargin)
    t = linspace(0, 2 * pi);

    x = c(1) + a*cos(t)*cos(theta) - b*sin(t)*sin(theta);
    y = c(2) + a*cos(t)*sin(theta) + b*sin(t)*cos(theta);

    [varargout{1:nargout}] = fill(x,y, varargin{:});
end