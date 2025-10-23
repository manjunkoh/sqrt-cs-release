function fig = plot_dispersion(mu, Sigma, samples, options)
    % given a mean (6 by 1) and covariance matrix (6 by 6), plot the samples and the ellipse for a lower triangular plot
    % samples is a 6 by N matrix

    arguments
        mu (6,1) double
        Sigma (6,6) double
        samples (6,:) double
        options.target_mu (6,1) double = NaN(6,1)
        options.target_Sigma (6,6) double = NaN(6,6)
        options.units = {'[LU]', '[VU]'}
    end
    target_mu = options.target_mu;
    target_Sigma = options.target_Sigma;
    units = options.units;

    % Sample statistics
    muMC = mean(samples, 2);
    SigmaMC = cov(samples');

    labels = {strcat('$x$, ', units(1)), strcat('$y$, ', units(1)), strcat('$z$, ', units(1)), strcat('$\dot{x}$, ', units(2)), strcat('$\dot{y}$, ', units(2)), strcat('$\dot{z}$, ', units(2))};
    lims = NaN(6,2);
    plot_scale = 7;
    for i = 1:6
        lims(i,:) = [mu(i) - plot_scale*sqrt(Sigma(i,i))  mu(i) + plot_scale*sqrt(Sigma(i,i))];
    end     

    fig = figure(WindowStyle='normal');
    set(fig, 'Position', [3 3 20 20]);
    t = tiledlayout(6,6);
    t.TileSpacing = 'none';
    t.Padding = 'compact'; % so that xlabel and exponents don't overlap
    mygreen = [62 150 81]/255;
   
    % TODO: make axis equal for position-position and velocity-velocity plots
    NumBins = 30;
    for i = 1:6
        for j = 1:6
            nexttile
            if i == j
                hold on
                histogram(samples(i,:), 'FaceColor', [0.5,0.5,0.5], 'NumBins', NumBins, 'DisplayName', 'Samples')
                normpdf = @(x) 1/sqrt(2*pi*Sigma(i,i)) * exp(-0.5*(x-mu(i)).^2/Sigma(i,i));
                normpdfMC = @(x) 1/sqrt(2*pi*SigmaMC(i,i)) * exp(-0.5*(x-muMC(i)).^2/SigmaMC(i,i));
                normpdfTarget = @(x) 1/sqrt(2*pi*target_Sigma(i,i)) * exp(-0.5*(x-target_mu(i)).^2/target_Sigma(i,i));
                x = linspace(lims(i,1), lims(i,2), 100);
                plot(x, normpdf(x)*length(samples(i,:))*(lims(i,2)-lims(i,1))/NumBins, 'r', 'DisplayName', 'Prediction Distribution/3$\sigma$ Ellipse')
                plot(x, normpdfMC(x)*length(samples(i,:))*(lims(i,2)-lims(i,1))/NumBins, 'b', 'DisplayName', 'Sample Distribution/3$\sigma$ Ellipse')
                if ~isnan(target_mu)
                    plot(x, normpdfTarget(x)*length(samples(i,:))*(lims(i,2)-lims(i,1))/NumBins, 'Color', mygreen, 'DisplayName', 'Target Distribution/3$\sigma$ Ellipse')
                end
                hold off
                xlim(lims(i,:))
                axis square
            elseif i > j
                hold on
                scatter(samples(j,:), samples(i,:), '.',  'MarkerEdgeColor', [0.5,0.5,0.5], 'MarkerFaceColor', [0.5,0.5,0.5],  DisplayName='Samples')
                plot3sigmaEllipse(mu([j,i]), Sigma([j,i],[j,i]), 'r', DisplayName='Prediction 3$\sigma$ Ellipse')
                plot3sigmaEllipse(muMC([j,i]), SigmaMC([j,i],[j,i]), 'b', DisplayName='Sample 3$\sigma$ Ellipse')
                if ~isnan(target_mu)
                    plot3sigmaEllipse(target_mu([j,i]), target_Sigma([j,i],[j,i]), Color=mygreen, DisplayName='Target 3$\sigma$ Ellipse')
                end
                hold off
                xlim(lims(j,:))
                ylim(lims(i,:))
                axis square
                % axis equal
            else
                axis off
            end
            if i == 6
                xlabel(labels{j}, 'FontSize', 20);
            else
                xticklabels({})
            end
            if j == 1 && i ~= 1
                ylabel(labels{i}, 'FontSize', 20);
            else
                yticklabels({})
            end
        end
    end
    
    lgd = legend('FontSize', 20);
    lgd.Layout.Tile = 4;

    figure_settings('medium');

end