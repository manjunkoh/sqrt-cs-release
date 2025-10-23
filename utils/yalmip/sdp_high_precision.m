function options = sdp_high_precision(options)
    if nargin < 1
        options = sdpsettings('verbose', 0, 'solver', 'mosek');
    end

    switch options.solver
        case 'mosek'
            options = mosek_high_precision(options);
        case 'sdpt3'
            options = sdpt3_high_precision(options);
        case 'sedumi'
            options = sedumi_high_precision(options);
        case 'copt'
            options = copt_high_precision(options);
        case 'sdpa_gmp'
            options = sdpa_gmp_high_precision(options);
        otherwise
            error('Unknown solver');
    end

    options = mosek_high_precision(options);
    options = sdpt3_high_precision(options);
    options = sedumi_high_precision(options);
    options = copt_high_precision(options);
    options = sdpa_gmp_high_precision(options);    

end


function options = mosek_high_precision(options)
% Mosek parameter documentation
% https://docs.mosek.com/latest/cxxfusion/parameters.html
% CVX precision effect on Mosek parameters
% https://groups.google.com/g/mosek/c/QMM_Bsx44AI
% https://ask.cvxr.com/t/pro-tip-mosek-and-cvx-precision/8364
% https://docs.mosek.com/latest/dotnetapi/solving-conic.html#adjusting-optimality-criteria

% options.mosek.MSK_DPAR_ANA_SOL_INFEAS_TOL = 1e-10; % Not in Mosek docs, but added in Yalmip
% options.mosek.MSK_DPAR_BASIS_TOL_S        = 1e-9; % default 1e-6,      min 1e-9
% options.mosek.MSK_DPAR_BASIS_TOL_X        = 1e-9; % default 1e-6,      min 1e-9

% options.mosek.MSK_DPAR_INTPNT_TOL_STEP_SIZE = 1e-12; % default 1e-6, min 0


% options.mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED  = 1e-10; % default 1e-8,  min 0
options.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-10; % default 1e-8,  min 0
options.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-10; % default 1e-8,  min 0
options.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-10; % default 1e-8,  min 0
% options.mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS  = 1e-12; % default 1e-12, min 0

% options.mosek.MSK_DPAR_INTPNT_TOL_PATH = 0.01; % default 1e-8, [0, 0.9999]

% options.mosek.MSK_DPAR_INTPNT_CO_TOL_NEAR_REL                  = 1; % default 1000

% presolve parameters
% options.mosek.MSK_DPAR_PRESOLVE_TOL_ABS_LINDEP                 = 1e-12; % default 1e-6, min 0
% options.mosek.MSK_DPAR_PRESOLVE_TOL_AIJ                        = 1e-12; % default 1e-12, min 1e-15
% options.mosek.MSK_DPAR_PRESOLVE_TOL_PRIMAL_INFEAS_PERTURBATION = 1e-12; % default 1e-6, min 0
% options.mosek.MSK_DPAR_PRESOLVE_TOL_REL_LINDEP                 = 1e-12; % default 1e-10, min 0
% options.mosek.MSK_DPAR_PRESOLVE_TOL_S                          = 1e-12; % default 1e-8, min 0
% options.mosek.MSK_DPAR_PRESOLVE_TOL_X                          = 1e-12; % default 1e-8, min 0

% maximum iterations for interior point method
options.mosek.MSK_IPAR_INTPNT_MAX_ITERATIONS = 2000; % default 400, max inf

end

function options = sdpt3_high_precision(options)

    options.sdpt3.gaptol = 1E-12;
    options.sdpt3.inftol = 1E-12;
    options.sdpt3.maxit = 400;

end

function options = sedumi_high_precision(options)

    options.sedumi.eps = 1E-12;
    options.sedumi.numtol = 1E-12;
end

function options = copt_high_precision(options)

    options.copt.FeasTol = 1E-12;
    options.copt.DualTol = 1E-12;
    options.copt.RelGap = 1E-12;
    options.copt.AbsGap = 1E-12;

end

function options = sdpa_gmp_high_precision(options)
    % sdpa-gmp can only be used on Linux/Unix
    % See for default parameters:
    % https://github.com/nakatamaho/sdpa-gmp/blob/master/param.sdpa
    % See for parameter documentation:
    % https://www.researchgate.net/publication/247456489_SDPA_SemiDefinite_Programming_Algorithm_User's_Manual_-_Version_600
    options.sdpa_gmp.maxIteration = 2000;
    options.sdpa_gmp.precision = 30;
    options.sdpa_gmp.epsilonDash = 1E-12;
    options.sdpa_gmp.epsilonStar = 1E-12;
    % options.sdpa_gmp.lowerBound = -1E-5;
    % options.sdpa_gmp.upperBound = 1E-5;
    % options.sdpa_gmp.lambdaStar = 1E4;
    % options.sdpa_gmp.betaBar = 0.3;
    % options.sdpa_gmp.gammaStar = 0.9;

end