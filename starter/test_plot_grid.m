L = 20;
grid = create_grid(L);
disp('Grid created');

scatter3(grid(:, 1), grid(:, 2), grid(:, 3), 10, grid(:, 3));
axis equal;
