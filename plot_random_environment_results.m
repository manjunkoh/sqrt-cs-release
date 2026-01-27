
trial = 1;
figure(Position=[0, 0, 15, 15])
sgtitle(sprintf('Trial %d', trial))

tiledlayout(2, 1);
nexttile;
title(sprintf('Full Covariance, Objective: %.3f', results.objective_full_covariance(trial)))
hold on
plot_problem(results.obstacle_centers{trial}, results.obstacle_radii{trial}, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);

if ~isempty(results.prob_full_covariance{trial}) && ~isempty(results.prob_full_covariance{trial}.mu)
	plot_solution(results.prob_full_covariance{trial}.mu, results.prob_full_covariance{trial}.P, state_risk);
end
plot(results.x_opt_deterministic{trial}(1,:), results.x_opt_deterministic{trial}(2,:), 'r.-', DisplayName='deterministic', LineWidth=1);


xlim([-1.4, 11])
ylim([-3, 2.0])

nexttile;
title(sprintf('Sqrt QR, Objective: %.3f', results.objective_qr(trial)))
hold on
plot_problem(results.obstacle_centers{trial}, results.obstacle_radii{trial}, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);

if ~isempty(results.prob_qr{trial}) && ~isempty(results.prob_qr{trial}.mu)
	plot_solution(results.prob_qr{trial}.mu, results.prob_qr{trial}.P, state_risk);
end
plot(results.x_opt_deterministic{trial}(1,:), results.x_opt_deterministic{trial}(2,:), 'r.-', DisplayName='deterministic', LineWidth=1);

% legend(legendUnq(), Location='northoutside', Orientation='horizontal', IconColumnWidth=15, FontSize=25)
xlim([-1.4, 11])
ylim([-3, 2.0])

save_filename = sprintf('./data/obstacle_planning_test_%s/trial_%d.pdf', current_time, trial);

