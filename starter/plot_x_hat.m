% Look at:
% 1. http://math.boisestate.edu/~wright/montestigliano/PlottingOnTheSphere.pdf
%    for info on how to plot on spheres in various ways
% 2. http://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot
%    for creating subplots with easily configurable spacing
% 3. https://www.mathworks.com/matlabcentral/newsreader/view_thread/305046
%    for interpolation on a sphere
% 4. http://www.mathworks.com/help/matlab/ref/griddata.html
%    Reference for griddata
% 5. https://www.mathworks.com/matlabcentral/newsreader/view_thread/59500
%    Issue with griddata not being able to interpolate outside the convex hull
% 6. http://www.mathworks.com/help/matlab/math/interpolating-scattered-data.html
%    Matlab help page on interpolation from scattered points

function plot_x_hat(x_hat, thetas, phis)

%l_max = 22;
%load(sprintf('grid_L%d.mat', l_max));
%load(fname);

%x_hat = x_hat';

%x_hat = zeros(size(x_hat));
%x_hat(7^2 + 7 + 1) = 1;

%v = real(mean(x_hat, 2));
v = x_hat;

%X = (r .* sin(thetas) .* cos(phis))';
%Y = (r .* sin(thetas) .* sin(phis))';
%Z = (r .* cos(thetas))';

res = 201;
t = linspace(0, pi, ceil(res/2));
p = linspace(-pi, pi, res);
[T, P] = meshgrid(t, p);

vq = griddata(thetas, phis, v, T, P, 'v4');

% Cortex
r = 8;
x = r .* sin(T) .* cos(P);
y = r .* sin(T) .* sin(P);
z = r.* cos(T);

%figure;
surf(x, y, z, vq);
%colorbar;
colormap jet;
caxis([min(vq(:)), max(vq(:))]);
shading interp;
daspect([1 1 1]);
hold on;

c = 0.9 * [1, 1, 1];

% Skull
r = 8.7;
x = r .* sin(T) .* cos(P);
y = r .* sin(T) .* sin(P);
z = r.* cos(T);
skull = surf(x, y, z);
set(skull, 'FaceColor', c, 'FaceAlpha', 0.1, 'EdgeColor', c, 'EdgeAlpha', 0.1);
shading interp;
hold on;

% Scalp
r = 9.2;
x = r .* sin(T) .* cos(P);
y = r .* sin(T) .* sin(P);
z = r.* cos(T);
scalp = surf(x, y, z);
set(scalp, 'FaceColor', c, 'FaceAlpha', 0.1, 'EdgeColor', c, 'EdgeAlpha', 0.1);
shading interp;

axis off;
axis tight;
vw = [70 25];
view(vw);
