function solve_flag = solve_subproblem_optimizer(obj)
    % Get the numerical values of the parameters for the optimizer
    lambda   = obj.this_iter.lambda;
    mu       = obj.this_iter.mu;
    w        = obj.this_iter.w;
    r        = obj.this_iter.r;
    ref_vars = obj.this_iter.ref_vars;

    changing_params = obj.scp_prob.get_changing_parameters(ref_vars);
    changing_params_cell = struct2cell(changing_params);

    all_params_double = [{lambda}, {mu}, {w}, {r}, changing_params_cell(:)'];
    all_params_double = all_params_double(~cellfun('isempty', all_params_double));
    for idx_wanted_vars = 1:length(obj.var_names)
        field = obj.var_names{idx_wanted_vars};
        if obj.scp_prob.is_used_in_constraint_update.(field)
            all_params_double{end+1} = ref_vars.(field);
        end
    end

    [wanted_vars, solve_flag, ~, ~, ~, diagnostics] = obj.optimizer(all_params_double);

    obj.this_iter.J0 = wanted_vars{1};
    obj.this_iter.L = wanted_vars{2};
    % For debugging, uncomment xi, xi_norm_squared, etc.
    idx_wanted_vars = 3;
    if obj.scp_prob.size_noncvx_eq > 0
        % xi = wanted_vars{i};
        idx_wanted_vars = idx_wanted_vars + 1;
        % xi_norm_squared = wanted_vars{i};
        idx_wanted_vars = idx_wanted_vars + 1;
    end
    if obj.scp_prob.size_noncvx_ineq > 0
        % zeta = wanted_vars{i};
        idx_wanted_vars = idx_wanted_vars + 1;
        % zeta_norm_squared = wanted_vars{i};
        idx_wanted_vars = idx_wanted_vars + 1;
    end
    obj.this_iter.vars = cell2struct(wanted_vars(idx_wanted_vars:end)', obj.var_names);
    obj.this_iter.DeltaZ = obj.getDeltaZ(obj.this_iter.vars);
    obj.this_iter.time_solver = diagnostics.solvertime;
    if isfield(diagnostics, 'yalmiptime')
        obj.this_iter.time_yalmip = diagnostics.yalmiptime;
    else
        obj.this_iter.time_yalmip = NaN;
    end
end