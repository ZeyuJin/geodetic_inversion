close all
clc
clear

addpath('/Users/zej011/work/Kang_tutorial/codes_utilities/matlab/igppsar');
addpath(genpath('/Users/zej011/work/zeyu/matlab/external_func'));


%% load fault surface trace
[xf1, yf1, ~] = read_data_from_file('~/Qinghai/geometry/curveFault/data/curved_cut.txt');
[xf2, yf2, ~] = read_data_from_file('~/Qinghai/geometry/curveFault/data/curved2.txt');
[xf3, yf3, ~] = read_data_from_file('~/Qinghai/geometry/curveFault/data/curved3.txt');

% linear interpolation of surface points
xq1 = linspace(min(xf1), max(xf1), length(xf1)*1e2)';
yq1 = interp1(xf1, yf1, xq1, 'spline');

xq2 = linspace(min(xf2), max(xf2), length(xf2)*1e2)';
yq2 = interp1(xf2, yf2, xq2, 'spline');

xq3 = linspace(min(xf3), max(xf3), length(xf3)*1e2)';
yq3 = interp1(xf3, yf3, xq3, 'spline');


%% incorporate seismicity catalog
[xs, ys, sdepth] = read_data_from_file('~/Qinghai/geometry/He_etal_data/aftershocks.txt');
% mask outliers
[mkx, mky, ~] = read_data_from_file('~/Qinghai/geometry/curveFault/data/aftershocks_mask.txt');

in = inpolygon(xs, ys, mkx, mky);  % within the polygon
xtmp = xs(in);
dtmp = sdepth(in);
ytmp = ys(in);

in2 = find(dtmp >= -25);  % depth < 20km
xs = xtmp(in2);
ys = ytmp(in2);
sdepth = dtmp(in2);
clear xtmp dtmp ytmp


%% fit the surface 
% apply larger weight for points at the surface
% The 'cubicinterp' does not need weight vector

% fit the surface segment by segment
% cut seismicity close to each segment
xtail = 0;   % 5km transition zone along east direction
leftEnd = xs >= min(xf1)-xtail;
rightEnd = xs <= max(xf1)+xtail;
depthCon1 = sdepth < -5;
depthCon2 = sdepth > -20;
xs1Cut = xs(leftEnd & rightEnd & depthCon1 & depthCon2);
ys1Cut = ys(leftEnd & rightEnd & depthCon1 & depthCon2);
depth1Cut = sdepth(leftEnd & rightEnd & depthCon1 & depthCon2);


% % weighted gridfit (also works)
% xOne = [xf1; xs1Cut];
% depthOne = [zeros(size(xf1)); depth1Cut];
% yOne = [yf1; ys1Cut];
% 
% % assign larger weights to surface data
% ratio = 1e2;
% w1 = ones(size(xOne));
% w1(1:length(xf1)) = ratio;


% bifucated tails
xtail = 5;
ind2 = find(xs >= min(xf2)-xtail & xs <= max(xf2)+xtail);
xs2Cut = xs(ind2);
ys2Cut = ys(ind2);
depth2Cut = sdepth(ind2);

% separate seismicity for each tail (mannually)
[l2x, l2y] = read_data_from_file('data/L2mask.txt');
l2In = inpolygon(xs2Cut, ys2Cut, l2x, l2y);
l1In = ~inpolygon(xs2Cut, ys2Cut, l2x, l2y);

xsCutSecond = xs2Cut(l2In);
ysCutSecond = ys2Cut(l2In);
depthCutSecond = depth2Cut(l2In);

xsCutFirst = xs2Cut(l1In);
ysCutFirst = ys2Cut(l1In);
depthCutFirst = depth2Cut(l1In);

% fitted points for each segment
xOne = [xq1; xs1Cut];
depthOne = [zeros(size(xq1)); depth1Cut];
yOne = [yq1; ys1Cut];

xTwo = [xq2; xsCutSecond];
depthTwo = [zeros(size(xq2)); depthCutSecond];
yTwo = [yq2; ysCutSecond];

xThree = [xq3; xsCutFirst];
depthThree = [zeros(size(xq3)); depthCutFirst];
yThree = [yq3; ysCutFirst]; 


%% create grid points
% mini = 0.75;  % in km, the minimum size along strike and depth
mini = 3;
maxD = -25;  % largest depth along the dip direction

% uniform depth distribution
nd = floor(abs(maxD)/mini);
% nd = 10;
d = linspace(maxD, 0, nd);

% first plane
nx1 = ceil((max(xf1) - min(xf1)) / mini);
% nx1 = 75;
x1 = linspace(min(xf1), max(xf1), nx1);
[x1m, d1m, y1m, x1p, d1p, y1p] = gridFitInterpolate(xOne, depthOne, yOne, x1, d);

% enforce a match between fault branches at the Y intersection
x1end = x1m(:,end);
d1end = d1m(:,end);
y1end = y1m(:,end);

xTwo = [x1end; xTwo];
depthTwo = [d1end; depthTwo];
yTwo = [y1end; yTwo];

xThree = [x1end; xThree];
depthThree = [d1end; depthThree];
yThree = [y1end; yThree];

% second plane
nx2 = ceil((max(xf2) - min(xf2)) / mini);
% nx2 = 15;
x2 = linspace(min(xf2), max(xf2), nx2);
[x2m, d2m, y2m, x2p, d2p, y2p] = gridFitInterpolate(xTwo, depthTwo, yTwo, x2, d);

% third plane 
nx3 = ceil((max(xf3) - min(xf3)) / mini);
% nx3 = 20;
x3 = linspace(min(xf3), max(xf3), nx3);
[x3m, d3m, y3m, x3p, d3p, y3p] = gridFitInterpolate(xThree, depthThree, yThree, x3, d);

% save the data
save('pointSourceGeometry/Curve1Mesh.mat','x1m','y1m','d1m');
save('pointSourceGeometry/Curve2Mesh.mat','x2m','y2m','d2m');
save('pointSourceGeometry/Curve3Mesh.mat','x3m','y3m','d3m');

save('pointSourceGeometry/Curve1Point.mat','x1p','y1p','d1p');
save('pointSourceGeometry/Curve2Point.mat','x2p','y2p','d2p');
save('pointSourceGeometry/Curve3Point.mat','x3p','y3p','d3p');

%% plot the surface
figure; hold on
plot3(xf1, zeros(size(xf1)), yf1, 'linewidth', 5, 'color', 'm');
plot3(xf2, zeros(size(xf2)), yf2, 'linewidth', 5, 'color', 'm');
plot3(xf3, zeros(size(xf3)), yf3, 'linewidth', 5, 'color', 'm');

% plot the seismicity
sz = 10;
scatter3(xs1Cut, depth1Cut, ys1Cut, sz, [0 0.5 0], 'filled');
scatter3(xs2Cut, depth2Cut, ys2Cut, sz, [0 0.5 0], 'filled');
scatter3(xTwo, depthTwo, yTwo, sz, [0 0.5 0], 'filled');
scatter3(xThree, depthThree, yThree, sz, [0 0.5 0], 'filled');

mesh( x1m, d1m, y1m, ...
    'LineStyle', '-', 'LineWidth', 2, 'EdgeColor', 'k', ...
    'FaceColor', 'none', 'FaceAlpha', 1 );

mesh( x2m, d2m, y2m, ...
    'LineStyle', '-', 'LineWidth', 2, 'EdgeColor', 'r', ...
    'FaceColor', 'none', 'FaceAlpha', 1 );

mesh( x3m, d3m, y3m, ...
    'LineStyle', '-', 'LineWidth', 2, 'EdgeColor', 'b', ...
    'FaceColor', 'none', 'FaceAlpha', 1 );

% % plot the point source distribution
% scatter3(x1p(:), d1p(:), y1p(:), sz, 'ko', 'filled');
% scatter3(x2p(:), d2p(:), y2p(:), sz, 'ro', 'filled');
% scatter3(x3p(:), d3p(:), y3p(:), sz, 'bo', 'filled');

set(gca, 'YDir','reverse');
axis equal
xlabel( 'East (km)' );
ylabel( 'Depth (km)' );
zlabel( 'North (km)' );
ylim([-25 0]);
xlim([min(xOne) 30]);
set(gca, 'fontsize', 20);
