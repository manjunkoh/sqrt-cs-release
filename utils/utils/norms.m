function cvx_optval = norms( x, p)

    % NORMS: A lightweight version of the NORMS function included with CVX.
if nargin < 2 || isempty( p )
    p = 2;
elseif ~isnumeric( p ) || numel( p ) ~= 1 || ~isreal( p )
    error( 'Second argument must be a real number.' );
elseif p < 1 || isnan( p )
    error( 'Second argument must be between 1 and +Inf, inclusive.' );
end
    
dim = 1;

switch p
    case 1
        cvx_optval = sum( abs( x ), dim );
    case 2
        cvx_optval = sqrt( sum( x .* x, dim ) );
    case Inf
        cvx_optval = max( abs( x ), [], dim );
    otherwise
        cvx_optval = sum( abs( x ) .^ p, dim ) .^ ( 1 / p );
end