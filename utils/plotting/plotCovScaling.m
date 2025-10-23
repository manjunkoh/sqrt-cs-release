function [max_magnitudes, min_magnitudes] = plotCovScaling(p, A,B,Q,m)

    % Q(:,:,k) * vars.Pest(:,:,k+1) * Q(:,:,k)'... 
    %         - Q(:,:,k) * A(:,:,k)*vars.Pest(:,:,k)* A(:,:,k)' *Q(:,:,k)'...
    %         - Q(:,:,k) * A(:,:,k)*vars.U(:,:,k)'* B(:,:,k)'*Q(:,:,k)'...
    %         - Q(:,:,k) * B(:,:,k)*vars.U(:,:,k)*A(:,:,k)'* Q(:,:,k)'...
    %         - Q(:,:,k) * B(:,:,k)*vars.Y(:,:,k)*B(:,:,k)'* Q(:,:,k)'...
    %         - Q(:,:,k) * m.KalGain(:,:,k) * m.Pymintil(:,:,k) * m.KalGain(:,:,k)'*Q(:,:,k)' == 0
    
    num_terms = 6;
    max_magnitudes = zeros(num_terms, p.Nseg);
    for k = 1:p.Nseg
        max_magnitudes(1,k) = max(abs(Q(:,:,k) * ones(6,6) * Q(:,:,k)'), [], 'all');
        max_magnitudes(2,k) = max(abs(Q(:,:,k) * A(:,:,k) * ones(6,6) * A(:,:,k)' * Q(:,:,k)'), [], 'all');
        max_magnitudes(3,k) = max(abs(Q(:,:,k) * A(:,:,k) * ones(6,3) * B(:,:,k)' * Q(:,:,k)'), [], 'all');
        max_magnitudes(4,k) = max(abs(Q(:,:,k) * B(:,:,k) * ones(3,6) * A(:,:,k)' * Q(:,:,k)'), [], 'all');
        max_magnitudes(5,k) = max(abs(Q(:,:,k) * B(:,:,k) * ones(3,3) * B(:,:,k)' * Q(:,:,k)') , [], 'all');
        max_magnitudes(6,k) = max(abs(Q(:,:,k) * m.KalGain(:,:,k) * m.Pymintil(:,:,k) * m.KalGain(:,:,k)' * Q(:,:,k)'), [], 'all');
    end

    min_magnitudes = zeros(num_terms, p.Nseg);
    for k = 1:p.Nseg
        min_magnitudes(1,k) = min(abs(Q(:,:,k) * ones(6,6) * Q(:,:,k)'), [], 'all');
        min_magnitudes(2,k) = min(abs(Q(:,:,k) * A(:,:,k) * ones(6,6) * A(:,:,k)' * Q(:,:,k)'), [], 'all');
        min_magnitudes(3,k) = min(abs(Q(:,:,k) * A(:,:,k) * ones(6,3) * B(:,:,k)' * Q(:,:,k)'), [], 'all');
        min_magnitudes(4,k) = min(abs(Q(:,:,k) * B(:,:,k) * ones(3,6) * A(:,:,k)' * Q(:,:,k)'), [], 'all');
        min_magnitudes(5,k) = min(abs(Q(:,:,k) * B(:,:,k) * ones(3,3) * B(:,:,k)' * Q(:,:,k)'), [], 'all');
        min_magnitudes(6,k) = min(abs(Q(:,:,k) * m.KalGain(:,:,k) * m.Pymintil(:,:,k) * m.KalGain(:,:,k)' * Q(:,:,k)'), [], 'all');
    end

    figure
    hold on
    C = colororder;

    for i = 1:num_terms
        plot(max_magnitudes(i,:)', Color=C(i,:), LineStyle='-', DisplayName=sprintf('Term %d', i))
        plot(min_magnitudes(i,:)', Color=C(i,:), LineStyle='--', HandleVisibility='off')
    end

    legend()
    xlabel('Segment')
    ylabel('Magnitude')
    title('Max and Min Magnitudes of Coefficient Matrices in Covariance Equation')
end

