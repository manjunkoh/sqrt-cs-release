function J = test_solve(objective, constraints, names)

    constraint_subset = constraint_list(constraints, names);
    ops = sdpsettings('solver', 'mosek', 'verbose', 0);
    sol = optimize(constraint_subset, objective, ops);

    if sol.problem == 0
        disp('YALMIP: Problem solved')
        J = value(objective);
    else
        disp(yalmiperror(sol.problem))
        J = NaN;
    end
end
