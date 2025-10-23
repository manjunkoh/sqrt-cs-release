function varargout = plotErrorEllipse(mu, cov, p, varargin)
    % https://www.xarg.org/2018/04/how-to-plot-a-covariance-error-ellipse/
    s = -2 * log(1 - p); %chi2inv(p,2)
    assert(all(eig(cov(1:2, 1:2))));
%     assert(isreal(eig(cov(1:2, 1:2))));
    [V,d] = eig(cov(1:2, 1:2));

    V1 = V(:,1);
    try
        theta = atan2(V1(2), V1(1)); % angle of the largest eigenvector with the x-axis
    catch
        warning('Small imaginary value in covariance matrix eigenvalue.');
        return
    end
    a = sqrt(s*d(1,1));
    b = sqrt(s*d(2,2));

%     plot(mu(1)+ semi_major*cos(t+alpha), mu(2)+ semi_minor*sin(t+alpha), 'Color', color)
    [varargout{1:nargout}] = plotEllipse(a,b,mu(1:2), theta, varargin{:});
    hold on
    return 
