function L_interp = interpolate_lower_triangular(L0, LN, N, method)
	%INTERPOLATE_LT Interpolates between two lower-triangular matrices L0 and LN
	assert(all(size(L0) == size(LN)), 'L0 and LN must have the same dimensions');
	assert(istril(L0), 'L0 must be lower-triangular');
	assert(istril(LN), 'LN must be lower-triangular');
	n = size(L0, 1);
	switch method
		case 'cholesky'
			L_interp = zeros(n, n, N);
			for i = 1:size(L0, 1)
				for j = 1:size(L0, 2)
					if i < j
						continue;
					end
					L_interp(i,j,:) = linspace(L0(i,j), LN(i,j), N);
				end
			end
		case 'log-cholesky'
			% interpolate the off-diagonal entries linearly
			L_interp = zeros(n, n, N);
			for i = 1:n
				for j = 1:n	
					if i <= j
						% Diagonal entries will be handled separately
						continue;
					end
					L_interp(i,j,:) = linspace(L0(i,j), LN(i,j), N);
				end
			end

			ts = linspace(0, 1, N);
			for k = 1:N
				t = ts(k);
				for i = 1:n
					L_interp(i,i,k) = L0(i,i) * exp(t * log(LN(i,i)/L0(i,i)));
				end
			end
		otherwise
			error('Unknown interpolation method: %s', method);
	end
end
				
			



