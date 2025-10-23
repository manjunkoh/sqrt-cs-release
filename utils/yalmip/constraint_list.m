function constraint_subset = constraint_list(constraints, names)

    constraint_subset = [];
    for i = 1:length(names)
        constraint_subset = [constraint_subset; constraints(names{i})];
    end
end
