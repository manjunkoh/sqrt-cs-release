function varargout = plot1sigmaEllipse(mu, cov, varargin)
%PLOT3SIGMAELLIPSE plots a 1 sigma ellipse

	assert(all(eig(cov(1:2, 1:2))));
	[V,d] = eig(cov(1:2, 1:2));

	V1 = V(:,1);
	try
		theta = atan2(V1(2), V1(1)); % angle of the largest eigenvector with the x-axis
	catch
		warning('Small imaginary value in covariance matrix eigenvalue.');
		return
	end

	s = 1;
	a = s * sqrt(d(1,1));
	b = s * sqrt(d(2,2));

	%     plot(mu(1)+ semi_major*cos(t+alpha), mu(2)+ semi_minor*sin(t+alpha), 'Color', color)
	[varargout{1:nargout}] = plotEllipse(a,b,mu(1:2), theta, varargin{:});

end