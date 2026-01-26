function [x_hist_all, u_hist_all] = simulate_samples(mu_0, Sigma_0, A_sys, B_sys, G_sys, K, mu, v, num_nodes, num_simulations)
    x_hist_all = zeros(size(mu, 1), num_nodes+1, num_simulations);
    u_hist_all = zeros(size(v, 1), num_nodes, num_simulations);
    for i = 1:num_simulations
        [x_hist, u_hist] = simulate_sample(mu_0, Sigma_0, A_sys, B_sys, G_sys, K, mu, v, num_nodes);
        x_hist_all(:,:,i) = x_hist;
        u_hist_all(:,:,i) = u_hist;
    end
end
