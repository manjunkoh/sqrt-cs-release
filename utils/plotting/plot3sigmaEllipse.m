function plot3sigmaEllipse(mu,cov, varargin)
%PLOT3SIGMAELLIPSE plots a 3 sigma ellipse

p = normcdf(3) - normcdf(-3);
plotErrorEllipse(mu,cov,p,varargin{:});

end