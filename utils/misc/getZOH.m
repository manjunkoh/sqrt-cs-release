function u = getZOH(t, t_his, u_his)
    
    % getZOH - get the value of the zero-order-hold control signal at time t
    %
    % Inputs:
    %    t - times at which to evaluate the control signal
    %    t_his - ZOH discretization time nodes (N+1 x 1)
    %    u_his - ZOH control signal values (n_u x N)
    % Outputs:
    %    u - control signals at time t

    u = zeros(size(u_his, 1), length(t));

    % For each time t, find the time segment right before t
    for i = 1:length(t)

        if t(i) < t_his(1) || t(i) > t_his(end)
            continue;
        end

        k = find(t_his <= t(i), 1, 'last');

        if k == length(t_his)
            u(:,i) = u_his(:,end);
        else
            u(:,i) = u_his(:,k);
        end
    end

end
    