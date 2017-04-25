%% Starter code for source localization on a real-brain model
clear;

load('headmodel.mat');
load('leadfield-real-26673-2050.mat');

[m, n] = size(L);
x = zeros(n, 1);
x(25) = 1;

y = L * x;

%% TSVD code goes here



%% Plot reconstruction, assuming solution is stored in x_hat

% Load skull surfaces
bnd = headmodel.vol.bnd;
wm = pial_wm_surfaces.wm;

% Plot the subsampled pial and white matter surfaces
f = figure;
%trisurf(wm.faces, wm.vertices(:, 1), wm.vertices(:, 2), wm.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.6, 0.4], 'FaceAlpha', 0.5);
trisurf(wm.faces, wm.vertices(:, 1), wm.vertices(:, 2), wm.vertices(:, 3), x_hat, 'FaceAlpha', 0.5);
hold on;

% Plot the brain, CSF, skull and scalp surfaces
trisurf(bnd(1).tri, bnd(1).pnt(:, 1), bnd(1).pnt(:, 2), bnd(1).pnt(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.6, 0.6, 0.6], 'FaceAlpha', 0.1);
trisurf(bnd(2).tri, bnd(2).pnt(:, 1), bnd(2).pnt(:, 2), bnd(2).pnt(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 0.1);
trisurf(bnd(3).tri, bnd(3).pnt(:, 1), bnd(3).pnt(:, 2), bnd(3).pnt(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.2, 0.2, 0.2], 'FaceAlpha', 0.1);
axis equal;
