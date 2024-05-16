close all
clc
clear

%% load data and models
% fault trace
[xf1, yf1, ~] = read_data_from_file('data/curved_cut.txt');
[xf2, yf2, ~] = read_data_from_file('data/curved2.txt');
[xf3, yf3, ~] = read_data_from_file('data/curved3.txt');

% seismicity
load('data/seism_mask.mat', 'xseism', 'yseism', 'dseism');

% curved point source locations
load('pointSourceGeometry/Curve1Mesh.mat','x1m','y1m','d1m');
load('pointSourceGeometry/Curve2Mesh.mat','x2m','y2m','d2m');
load('pointSourceGeometry/Curve3Mesh.mat','x3m','y3m','d3m');

% load('pointSourceGeometry/Curve1Point.mat','x1p','y1p','d1p');
% load('pointSourceGeometry/Curve2Point.mat','x2p','y2p','d2p');
% load('pointSourceGeometry/Curve3Point.mat','x3p','y3p','d3p');

%% plot
figure; hold on
% plot the surface trace
plot3(xf1, zeros(size(xf1)), yf1, 'linewidth', 5, 'color', 'm');
plot3(xf2, zeros(size(xf2)), yf2, 'linewidth', 5, 'color', 'm');
plot3(xf3, zeros(size(xf3)), yf3, 'linewidth', 5, 'color', 'm');

% plot the seismicity
sz = 25;
scatter3(xseism, dseism, yseism, sz, [0 0.5 0], 'filled');

% plot the mesh or point locations?
mesh( x1m, d1m, y1m, ...
    'LineStyle', '-', 'LineWidth', 2, 'EdgeColor', 'k', ...
    'FaceColor', 'none', 'FaceAlpha', 1 );

mesh( x2m, d2m, y2m, ...
    'LineStyle', '-', 'LineWidth', 2, 'EdgeColor', 'r', ...
    'FaceColor', 'none', 'FaceAlpha', 1 );

mesh( x3m, d3m, y3m, ...
    'LineStyle', '-', 'LineWidth', 2, 'EdgeColor', 'b', ...
    'FaceColor', 'none', 'FaceAlpha', 1 );

% scatter3(x1p(:), d1p(:), y1p(:), sz, 'ko', 'filled');
% scatter3(x2p(:), d2p(:), y2p(:), sz, 'ro', 'filled');
% scatter3(x3p(:), d3p(:), y3p(:), sz, 'bo', 'filled');

set(gca, 'YDir','reverse');
axis equal
xlabel( 'East (km)' );
ylabel( 'Depth (km)' );
zlabel( 'North (km)' );
ylim([-25 0]);
xlim([min(xf1) max(xf3)]);
set(gca, 'fontsize', 20);
set(gca, 'Color', 'None');

%% make movie
% % draw the mp4 movie
% videoFilename = fullfile(pwd,'pointSource_movie.mp4');
% vidfile = VideoWriter(videoFilename,'MPEG-4');
% vidfile.Quality = 100;
% %vidfile.Duration = 30;
% vidfile.FrameRate = 10;
% open(vidfile);

% i=0;
% inv = 2;
% for j=-80:280
%     view(i,j)
%     drawnow
% %     print('-djpeg90','-r600','mod_3d')
% 
%     frame = getframe(gcf);
%     writeVideo(vidfile, frame);
% end
% 
% close(vidfile)
