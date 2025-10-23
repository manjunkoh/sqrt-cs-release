function vars = assignValues(vars)
    % Assign values to sdpvar objects after solving the optimization problem
    
    % Get field names of the structure
    fields = fieldnames(vars);

    % Loop through the fields
    for i = 1:length(fields)
        % Get the field name
        field = fields{i};
        
        vars.(field) = value(vars.(field));
    end
end