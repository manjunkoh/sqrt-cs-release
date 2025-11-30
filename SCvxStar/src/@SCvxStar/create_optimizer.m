function opt = create_optimizer(obj)
    fprintf('Creating optimizer object...')
    time_optimizer_start = tic;
    sdp_ref_vars = obj.scp_prob.define_vars();

    w = sdpvar(1);
    constraints = [
        obj.scp_prob.sdp_convex_eq
        obj.scp_prob.sdp_convex_ineq
    ];

    xi = obj.scp_prob.slack_noncvx_eq;
    if obj.scp_prob.size_noncvx_eq > 0
        lambda = sdpvar(obj.scp_prob.size_noncvx_eq, 1);
        xi_norm_squared = sdpvar(1);

        % See Problem 4.26 in Boyd's book
        constraints = [constraints
            % norm([2 * xi; xi_norm_squared - 1]) <= xi_norm_squared + 1
            cone([2 * xi; xi_norm_squared - 1], xi_norm_squared + 1)
            xi_norm_squared >= 0
        ];
    else
        lambda = [];
        xi_norm_squared = 0;
    end

    zeta = obj.scp_prob.slack_noncvx_ineq;
    if obj.scp_prob.size_noncvx_ineq > 0
        mu = sdpvar(obj.scp_prob.size_noncvx_ineq, 1);
        zeta_norm_squared = sdpvar(1);

        constraints = [constraints
            % norm([2 * zeta; zeta_norm_squared - 1]) <= zeta_norm_squared + 1
            cone([2 * zeta; zeta_norm_squared - 1], zeta_norm_squared + 1)
            zeta >= 0
        ];
    else
        mu = [];
        zeta_norm_squared = 0;
    end

    J0 = obj.scp_prob.sdp_objective;
    J_aug = J0 + obj.penalty_function_optimizer(xi, zeta, lambda, mu, w, xi_norm_squared, zeta_norm_squared);

    r = sdpvar(1);

    trust_region_constraint = [obj.get_trust_region_constraint_lhs(obj.scp_prob.sdp_vars, sdp_ref_vars, r) <= 0];
    
    constraints = [constraints; trust_region_constraint];

    constraints = [constraints; obj.scp_prob.get_changing_constraints(obj.scp_prob.sdp_vars, sdp_ref_vars)];
    changing_params_cell = struct2cell(obj.scp_prob.changing_params);
    all_params = [{lambda}, {mu}, {w}, {r}, changing_params_cell(:)'];
    all_params = all_params(~cellfun('isempty', all_params));

    for i = 1:length(obj.var_names)
        field = obj.var_names{i};
        if obj.scp_prob.is_used_in_constraint_update.(field)
            all_params{end+1} = sdp_ref_vars.(field);
        end
    end

    vars_cell = struct2cell(obj.scp_prob.sdp_vars);
    wanted_vars = [{J0}, {J_aug}];
    if obj.scp_prob.size_noncvx_eq > 0
        wanted_vars = [wanted_vars, {xi}, {xi_norm_squared}];
    end
    if obj.scp_prob.size_noncvx_ineq > 0
        wanted_vars = [wanted_vars, {zeta}, {zeta_norm_squared}];
    end
    wanted_vars = [wanted_vars, vars_cell(:)'];
    wanted_vars = wanted_vars(~cellfun('isempty', wanted_vars));
    wanted_vars = wanted_vars(~cellfun(@isnumeric, wanted_vars));


    % If problem is parameterized, add the parameters to the optimizer
    % For example, if the problem has initial and final states as parameters
    % TODO: case where changing constraints are parameterized
    [problem_params, parameterized_constraints] = obj.scp_prob.get_parameterized_constraints(obj.scp_prob.sdp_vars);

    if isempty(problem_params)
        opt = optimizer(constraints, J_aug, obj.constParams.yalmip_options, all_params, wanted_vars);
    else
        all_params = [all_params, struct2cell(problem_params)'];
        constraints = [constraints; parameterized_constraints];
        num_problem_params = length(fieldnames(problem_params));
        num_other_params = length(all_params) - num_problem_params;
        opt = optimizer(constraints, J_aug, obj.constParams.yalmip_options, all_params, wanted_vars);
        problem_params_cell = struct2cell(obj.scp_prob.p_double);
        opt = opt(cell(1, num_other_params), problem_params_cell{:});
    end
    fprintf('Done! (took %f seconds)\n', toc(time_optimizer_start));

end