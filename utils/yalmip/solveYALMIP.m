function [diagnostics, solve_flag] = solveYALMIP(constraints, objective, yalmip_options)
    arguments
        constraints
        objective
        yalmip_options = sdpsettings();
    end
    diagnostics = optimize(constraints, objective, yalmip_options);
    solve_flag = diagnostics.problem;
    switch solve_flag
        case 0 % Solved

        case 1 % Infeasible. resolve problem without objective 
            fprintf('not solved: %s\n', diagnostics.info)
            feasibility_problem = optimize(constraints, [], yalmip_options);
            if feasibility_problem.problem == 0
                fprintf('Feasibility problem solved. There may be a problem in the objective. \n')                        
            else
                fprintf('Feasibility problem not solved. The problem is infeasible. \n')
            end
        case 4 % Numerical problems
            fprintf('Numerical problems: MOSEK Error Code %i\n', diagnostics.solveroutput.r)
        otherwise
            fprintf('not solved: %s\n', diagnostics.info)
    end
    if yalmip_options.savesolveroutput && ~isempty(diagnostics.solveroutput.res.rmsg)
       fprintf(diagnostics.solveroutput.res.rmsg)
    end

end