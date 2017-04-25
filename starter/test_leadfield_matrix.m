dipole_L = 25;
dipole_grid = create_grid(dipole_L);

sensor_L = 15;
sensor_grid = create_grid(sensor_L);

A = leadfield_matrix(dipole_grid, sensor_grid);

x = zeros(dipole_L^2, 1);
x(517) = 1;
y = A * x;

figure;
scatter3(dipole_grid(:, 1), dipole_grid(:, 2), dipole_grid(:, 3), 10, x);
figure;
scatter3(sensor_grid(:, 1), sensor_grid(:, 2), sensor_grid(:, 3), 10, y);
