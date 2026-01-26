function plot_simulations(x_hist_all)
    for i = 1:size(x_hist_all, 3)
        x_hist = x_hist_all(:,:,i);
        plot(x_hist(1,:), x_hist(2,:), Color=[0, 0, 0, 0.5], LineWidth=0.5);
    end
end

