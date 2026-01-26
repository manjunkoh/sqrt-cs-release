function plot_solution(mu, P, state_risk)
    for k = 1:size(mu, 2)
        fill3sigmaEllipse(mu(:,k), P(:,:,k), '', FaceColor='#0082B2', FaceAlpha=0.5, EdgeColor='none', DisplayName="$3 \sigma$ ellipse");
    end
    plot(mu(1,:), mu(2,:), 'k.-', DisplayName='mean', LineWidth=1);
end
