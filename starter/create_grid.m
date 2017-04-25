function [grid, thetas, phis] = create_grid(L)
%% Function to create an "optimal-sampling" grid with L^2 points

thetas = zeros(L, 1);
for i = 1:L
	if mod(i, 2) == 0
		theta = pi * (i - 1) / (2 * L - 1);
	else
		theta = pi * (2 * L - i) / (2 * L - 1);
	end
	thetas(i) = theta;
end

phis = zeros(L^2, 1);
for i = 0:L-1
	for j = -i:i
		phis(i^2 + i + j + 1) = 2 * pi * j / (2*i + 1);
	end
end

% Compute cartesian coordinates
thetas = l_to_lm(thetas);
x = sin(thetas) .* cos(phis);
y = sin(thetas) .* sin(phis);
z = cos(thetas);

% Concatenate to create grid
grid = [x, y, z];
