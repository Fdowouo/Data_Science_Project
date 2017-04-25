function A = leadfield_matrix(dipole_grid, sensor_grid)

% Compute cos(theta), for all angles theta between every dipole and sensor pair
dot_prod_matrix = sensor_grid * dipole_grid';
costheta_matrix = bsxfun(@rdivide, dot_prod_matrix, sqrt(sum(dipole_grid.^2, 2))');
costheta_matrix = bsxfun(@rdivide, costheta_matrix, sqrt(sum(sensor_grid.^2, 2)));

% Correct errors caused by limited floating-point precision, which result
% in |costheta| > 1 for some elements.
costheta_matrix(costheta_matrix > 1) = 1;
costheta_matrix(costheta_matrix < -1) = -1;
[num_sensors, num_dipoles] = size(costheta_matrix);

phi = zeros(numel(costheta_matrix), 1);

lmax = 30;
Cl = zdipole_leadfield(7.5);  % Compute Clm for dipole depth of 0.5mm.
Cl = Cl(1:lmax);
Clm = zeros(1, lmax^2);
l = 0:lmax-1;
Clm(l.^2+l+1) = Cl;

Ylm = compute_Ylm(lmax, costheta_matrix(:), phi);

A = reshape(Clm * Ylm, num_sensors, num_dipoles);
