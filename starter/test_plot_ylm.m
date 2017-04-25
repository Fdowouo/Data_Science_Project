L = 20;
[grid, thetas, phis] = create_grid(L);

ylm = compute_Ylm(10, cos(thetas), phis);

plot_lm = 6;

figure;
scatter3(grid(:, 1), grid(:, 2), grid(:, 3), 10, real(ylm(plot_lm, :)));
colorbar;
axis equal;

figure;
scatter3(grid(:, 1), grid(:, 2), grid(:, 3), 10, imag(ylm(plot_lm, :)));
colorbar;
axis equal;
