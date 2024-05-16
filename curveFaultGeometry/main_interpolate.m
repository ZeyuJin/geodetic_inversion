close all
clc
clear

addpath(genpath('~/work/zeyu/matlab'));
addpath('~/work/Kang_tutorial/edcmp_matlab');

% load point sources geometry
load('pointSourceGeometry/Curve1Mesh.mat','x1m','y1m','d1m');
load('pointSourceGeometry/Curve2Mesh.mat','x2m','y2m','d2m');
load('pointSourceGeometry/Curve3Mesh.mat','x3m','y3m','d3m');

load('pointSourceGeometry/Curve1Point.mat','x1p','y1p','d1p');
load('pointSourceGeometry/Curve2Point.mat','x2p','y2p','d2p');
load('pointSourceGeometry/Curve3Point.mat','x3p','y3p','d3p');

% count the point source ID using a row-major order
meshG1 = [reshape(x1m',[],1), reshape(d1m',[],1), reshape(y1m',[],1)];
pointG1 = [reshape(x1p',[],1), reshape(d1p',[],1), reshape(y1p',[],1)];

meshG2 = [reshape(x2m',[],1), reshape(d2m',[],1), reshape(y2m',[],1)];
pointG2 = [reshape(x2p',[],1), reshape(d2p',[],1), reshape(y2p',[],1)];

meshG3 = [reshape(x3m',[],1), reshape(d3m',[],1), reshape(y3m',[],1)];
pointG3 = [reshape(x3p',[],1), reshape(d3p',[],1), reshape(y3p',[],1)];

% Original coarse tessellation
% biasL = 1.3;
% biasW = 1.2;

% geometric progression parameters
biasL = 1.4;
biasW = 1.2;
W = 25;

% The previous method that uses local points has some errors 
% use cross product to compute the exact strike and dip angles!

% % compute the area and interpolated on scattered points
% % original grid points from the bottom to the surfaces
% [strk1, dip1, ~] = compute_geometry_pointSource(x1m, y1m, d1m);
% [strk2, dip2, ~] = compute_geometry_pointSource(x2m, y2m, d2m);
% [strk3, dip3, ~] = compute_geometry_pointSource(x3m, y3m, d3m);
% 
% F1_strk = scatteredInterpolant(x1p(:), y1p(:), d1p(:), strk1(:));
% F2_strk = scatteredInterpolant(x2p(:), y2p(:), d2p(:), strk2(:));
% F3_strk = scatteredInterpolant(x3p(:), y3p(:), d3p(:), strk3(:));
% 
% F1_dip = scatteredInterpolant(x1p(:), y1p(:), d1p(:), dip1(:));
% F2_dip = scatteredInterpolant(x2p(:), y2p(:), d2p(:), dip2(:));
% F3_dip = scatteredInterpolant(x3p(:), y3p(:), d3p(:), dip3(:));


% interpolate uniform grids using depth-dependent progression
% ratio = 2;
ratio = 3;
[Vx1, Vy1, Vd1] = depth_dependent_point(meshG1, pointG1, biasL, biasW, W, ratio);
[Vx2, Vy2, Vd2] = depth_dependent_point(meshG2, pointG2, biasL, biasW, W, ratio);
[Vx3, Vy3, Vd3] = depth_dependent_point(meshG3, pointG3, biasL, biasW, W, ratio);


% build triangulation
% the row number of DT.ConnectivityList serves as the triangle IDs.
% save the info of bounding vertex
DT1 = delaunayTriangulation(Vx1, Vd1);
DT2 = delaunayTriangulation(Vx2, Vd2);
DT3 = delaunayTriangulation(Vx3, Vd3);

% build point source geometry using triangle vertices
miniScale = 1.5;  % minimum distance span between point sources

ntr1 = size(DT1.ConnectivityList, 1);
ntr2 = size(DT2.ConnectivityList, 1);
ntr3 = size(DT3.ConnectivityList, 1);
npara = ntr1+ntr2+ntr3;


% locate point sources into perspective triangles
% make sure that length(unique(ID1)) == size(DT1.ConnectivityList, 1)
% ID.shape = [number of point sources, 1]
% length(unique(ID)) = # of triangles
% neighbors.shape = [# of triangles, 3];
% NaN values in the neighbors for edge triangles

% Points (P1, P2, ...) are in row-major order
P1 = pointG1(:,1:2);
ID1 = pointLocation(DT1,P1);       % which triangle (row #) includes the data P
n1 = neighbors(DT1, (1:ntr1)');  % neighboring triangle ID (build smoothness matrix)

P2 = pointG2(:,1:2);
ID2 = pointLocation(DT2,P2);
n2 = neighbors(DT2, (1:ntr2)');

P3 = pointG3(:,1:2);
ID3 = pointLocation(DT3,P3);
n3 = neighbors(DT3, (1:ntr3)');

% % save the triangle info
% save('pointSourceGeometry/triangle_fault1.mat', 'ID1', 'n1', 'DT1', 'Vy1');
% save('pointSourceGeometry/triangle_fault2.mat', 'ID2', 'n2', 'DT2', 'Vy2');
% save('pointSourceGeometry/triangle_fault3.mat', 'ID3', 'n3', 'DT3', 'Vy3');


% plot the triangle
figure; hold on
trimesh(DT1.ConnectivityList, Vx1, Vd1, Vy1, 'linewidth', 2, 'EdgeColor', 'k', 'FaceColor', 'none');
trimesh(DT2.ConnectivityList, Vx2, Vd2, Vy2, 'linewidth', 2, 'EdgeColor', 'b', 'FaceColor', 'none');
trimesh(DT3.ConnectivityList, Vx3, Vd3, Vy3, 'linewidth', 2, 'EdgeColor', 'r', 'FaceColor', 'none');

% sz = 10;
% scatter3(x1p(:), d1p(:), y1p(:), sz, 'ko', 'filled'); 

% scatter3(v1(1), v1(3), v1(2), sz, 'ko', 'filled');
% scatter3(v2(1), v2(3), v2(2), sz, 'ko', 'filled');
% scatter3(v3(1), v3(3), v3(2), sz, 'ko', 'filled');
% trimesh(DT1.ConnectivityList(i,:), Vx1, Vd1, Vy1, 'linewidth', 2, 'EdgeColor', 'k', 'FaceColor', 'none');
% scatter3(tmpX, tmpZ, tmpY, sz-50, 'mo', 'filled');

% test smoothness
% trimesh(DT1.ConnectivityList(tri,:), Vx1, Vd1, Vy1, 'linewidth', 2, 'EdgeColor', 'k', 'FaceColor', 'none');
% trimesh(DT2.ConnectivityList(tr2-ntr1,:), Vx2, Vd2, Vy2, 'linewidth', 2, 'EdgeColor', 'k', 'FaceColor', 'none');
% trimesh(DT3.ConnectivityList(tr3,:), Vx3, Vd3, Vy3, 'linewidth', 2, 'EdgeColor', 'k', 'FaceColor', 'none');

% % test first triangle
% triID = 1;
% pId = find(ID1 == triID);
% 
% VertexID = DT1.ConnectivityList(triID,:);
% % trimesh(DT1.ConnectivityList(triID,:), Vx1, Vd1, Vy1, 'linewidth', 2, 'EdgeColor', 'r', 'FaceColor', 'none');
% % trimesh(DT1.ConnectivityList(ID1(1),:), Vx1, Vd1, Vy1, 'linewidth', 2, 'EdgeColor', 'g', 'FaceColor', 'none');
% 
% trimesh(DT1.ConnectivityList(138,:), Vx1, Vd1, Vy1, 'linewidth', 2, 'EdgeColor', 'k', 'FaceColor', 'none');
% trimesh(DT1.ConnectivityList(128,:), Vx1, Vd1, Vy1, 'linewidth', 2, 'EdgeColor', 'k', 'FaceColor', 'none');
% trimesh(DT1.ConnectivityList(15,:), Vx1, Vd1, Vy1, 'linewidth', 2, 'EdgeColor', 'k', 'FaceColor', 'none');
% trimesh(DT1.ConnectivityList(4,:), Vx1, Vd1, Vy1, 'linewidth', 2, 'EdgeColor', 'r', 'FaceColor', 'none');

% plot the fault trace
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


% % plot
% sz = 25;
% scatter3(x1p(:), d1p(:), y1p(:), sz-15, 'ko', 'filled'); 
% scatter3(P1(1,1), P1(1,2), pointG1(1,3), sz, 'go', 'filled');
% 
% scatter3(pointG1(pId,1), pointG1(pId,2), pointG1(pId,3), sz-15, 'ro', 'filled'); 
% % scatter3(x2p(:), d2p(:), y2p(:), sz, 'ko', 'filled');
% % scatter3(x3p(:), d3p(:), y3p(:), sz, 'ko', 'filled');
% % scatter3(Vx, Vd, Vy, sz, 'md', 'filled');
% sz = 5;
% scatter3(px, pz, py, 'ko', 'filled');

set(gca, 'YDir','reverse');
axis equal
ylim([-25 0]);
xlim([min(xf1) max(xf3)]);
view([0 45]);
set(gca, 'fontsize', 20);
set(gca, 'Color', 'None');
