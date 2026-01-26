function [x_hist, u_hist] = simulate_sample(mu_0, Sigma_0, A_sys, B_sys, G_sys, K, mu, v, num_nodes)
    x_hist = zeros(size(mu_0, 1), num_nodes+1);
    u_hist = zeros(size(v, 1), num_nodes);
    x_k = mu_0 + chol(Sigma_0, 'lower') * randn(size(mu_0));
    x_hist(:,1) = x_k;
    for k = 1:num_nodes
        u_k = v(:,k) + K(:,:,k) * (x_k - mu(:,k));
        x_k = A_sys(:,:,k) * x_k + B_sys(:,:,k) * u_k + G_sys(:,:,k) * randn(size(G_sys, 2), 1);

        x_hist(:,k+1) = x_k;
        u_hist(:,k) = u_k;
    end
end
