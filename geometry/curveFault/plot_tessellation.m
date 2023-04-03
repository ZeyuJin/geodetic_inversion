% close all
clc
clear

% load triangle info
load('pointSourceGeometry/triangle_fault1_original.mat', 'ID1', 'n1', 'DT1', 'Vy1');
load('pointSourceGeometry/triangle_fault2_original.mat', 'ID2', 'n2', 'DT2', 'Vy2');
load('pointSourceGeometry/triangle_fault3_original.mat', 'ID3', 'n3', 'DT3', 'Vy3');

Vx1 = DT1.Points(:,1);  Vd1 = DT1.Points(:,2);
Vx2 = DT2.Points(:,1);  Vd2 = DT2.Points(:,2);
Vx3 = DT3.Points(:,1);  Vd3 = DT3.Points(:,2);

% plot the triangle
figure; hold on
trimesh(DT1.ConnectivityList, Vx1, Vd1, Vy1, 'linewidth', 2, 'EdgeColor', 'k', 'FaceColor', 'none');
trimesh(DT2.ConnectivityList, Vx2, Vd2, Vy2, 'linewidth', 2, 'EdgeColor', 'b', 'FaceColor', 'none');
trimesh(DT3.ConnectivityList, Vx3, Vd3, Vy3, 'linewidth', 2, 'EdgeColor', 'r', 'FaceColor', 'none');

% % plot the fault trace
[xf1, yf1, ~] = read_data_from_file('~/Qinghai/geometry/curveFault/data/curved_cut.txt');
[xf2, yf2, ~] = read_data_from_file('~/Qinghai/geometry/curveFault/data/curved2.txt');
[xf3, yf3, ~] = read_data_from_file('~/Qinghai/geometry/curveFault/data/curved3.txt');

plot3(xf1, zeros(size(xf1)), yf1, 'linewidth', 5, 'color', 'm');
plot3(xf2, zeros(size(xf2)), yf2, 'linewidth', 5, 'color', 'm');
plot3(xf3, zeros(size(xf3)), yf3, 'linewidth', 5, 'color', 'm');

% plot the aftershocks
[xs, ys, sdepth] = read_data_from_file('~/Qinghai/geometry/He_etal_data/aftershocks.txt');
[mkx, mky, ~] = read_data_from_file('~/Qinghai/geometry/curveFault/data/aftershocks_mask.txt');

in = inpolygon(xs, ys, mkx, mky);  % within the polygon
xs = xs(in);
sdepth = sdepth(in);
ys = ys(in);

sz = 10;
scatter3(xs, sdepth, ys, sz, [0 0.5 0], 'filled');

% plot settings
set(gca, 'YDir','reverse');
ylabel('Depth (km)');
axis equal
xlabel('Easting (km)');
zlabel('Northing (km)');
ylim([-25 0]);
xlim([min(xf1) max(xf3)]);
view([0 45]);
set(gca, 'fontsize', 20);
set(gca, 'Color', 'None');
view([0 0]);
