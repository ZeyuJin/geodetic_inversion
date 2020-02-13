close all
clc
clear

addpath('/Users/zej011/work/Kang_tutorial/codes_utilities/matlab/igppsar');
addpath('/Users/zej011/work/Kang_tutorial/candis');
% addpath('/Users/zej011/Documents/MATLAB/matlab_functions/subaxis');

gps_file = '/Users/zej011/Ridgecrest/data_resample/GPS/total_fixed.txt';
this_track = '/Users/zej011/Ridgecrest/data_resample/ASC64/branch_cut';
insar_file = 'ASC64_0.2_cut.grd';
topo_file = 'samp_topo.grd';
east_angle_file = 'look_e.grd';
north_angle_file = 'look_n.grd';
up_angle_file = 'look_u.grd';

wavelength = 0.0554658*100; % C-band wavelength

[lonp,latp,ph] = grdread2([this_track,'/',insar_file]);
[lonp,latp,topo] = grdread2([this_track,'/',topo_file]);
topo = topo / 1000;
[lonp,latp,ze] = grdread2([this_track,'/',east_angle_file]);
[lonp,latp,zn] = grdread2([this_track,'/',north_angle_file]);
[lonp,latp,zu] = grdread2([this_track,'/',up_angle_file]);
[LON_INSAR,LAT_INSAR] = meshgrid(lonp,latp);

% downsample InSAR data for fast speed interpolation
dec = 2;
lonp = lonp(1:dec:end);
latp = latp(1:dec:end);
ph = ph(1:dec:end,1:dec:end);
topo = topo(1:dec:end,1:dec:end);
ze = ze(1:dec:end,1:dec:end);
zn = zn(1:dec:end,1:dec:end);
zu = zu(1:dec:end,1:dec:end);
LON_INSAR = LON_INSAR(1:dec:end,1:dec:end);
LAT_INSAR = LAT_INSAR(1:dec:end,1:dec:end);

gps_data = load(gps_file);
long = gps_data(:,1);
latg = gps_data(:,2);
de = gps_data(:,3);
dn = gps_data(:,4);
du = gps_data(:,5);

% find GPS within the bounds of interferogram
lonmin = min(lonp); lonmax = max(lonp);
latmin = min(latp); latmax = max(latp);
xv = [lonmin,lonmin,lonmax,lonmax,lonmin];
yv = [latmin,latmax,latmax,latmin,latmin];
gps_in = inpolygon(long,latg,xv,yv);
long_in = long(gps_in);
latg_in = latg(gps_in);
de = de(gps_in);
dn = dn(gps_in);
du = du(gps_in);

% Interpolated looking angles
Fe = scatteredInterpolant(double(LON_INSAR(:)),double(LAT_INSAR(:)),double(ze(:)));
gps_ze = Fe(double(long_in),double(latg_in));
Fn = scatteredInterpolant(double(LON_INSAR(:)),double(LAT_INSAR(:)),double(zn(:)));
gps_zn = Fn(double(long_in),double(latg_in));
Fu = scatteredInterpolant(double(LON_INSAR(:)),double(LAT_INSAR(:)),double(zu(:)));
gps_zu = Fu(double(long_in),double(latg_in));

% Interpolated phase and topo
Fph = scatteredInterpolant(double(LON_INSAR(:)),double(LAT_INSAR(:)),double(ph(:)));
unw = Fph(double(long_in),double(latg_in));
Ftopo = scatteredInterpolant(double(LON_INSAR(:)),double(LAT_INSAR(:)),double(topo(:)));
gps_topo = Ftopo(double(long_in),double(latg_in));

los_gps = (de .* gps_ze + dn .* gps_zn + du .* gps_zu) * 100;
los_insar = -1 * wavelength * unw / 4 / pi;
% los_insar = -1 * unw * 100;  % only for rng offsets
res = los_insar - los_gps;

% read the file with original resolution 
[lonp,latp,ph] = grdread2([this_track,'/',insar_file]);
[lonp,latp,topo] = grdread2([this_track,'/',topo_file]);
topo = topo / 1000;
[LON_INSAR,LAT_INSAR] = meshgrid(lonp,latp);

% fit a plane to the residual and subtract it from the interferogram
pfit = fit_ramp_topo(double(res),long_in,latg_in,double(gps_topo));
save([this_track,'/ramp_fit.txt'],'pfit','-ascii');        % for range offset use

los_ramp = pfit(1) .* LON_INSAR + pfit(2) .* LAT_INSAR + pfit(3) .* topo + pfit(4);
dis_insar = -1 * wavelength * ph / 4 / pi;
% dis_insar = -1 * ph * 100;
los_clean = dis_insar - los_ramp;

% plot boundary constrained by GPS sites
lonmin = min(lonp); lonmax = max(lonp);
latmin = min(latp); latmax = max(latp);

figure; 
sz = 30;
subplot('Position',[0.05 0.55 0.4 0.4]); hold on
pcolor(lonp,latp,dis_insar);
scatter(long,latg,sz,'k');
shading flat
colormap jet
colorbar
title('Raw InSAR Data');
axis([lonmin lonmax latmin latmax]);
set(gca,'Fontsize',20);

subplot('Position',[0.55 0.55 0.4 0.4]); hold on
pcolor(lonp,latp,los_clean);
scatter(long,latg,sz,'k');
shading flat
colormap jet
colorbar
title('Detrended InSAR Data');
axis([lonmin lonmax latmin latmax]);
set(gca,'Fontsize',20);

subplot('Position',[0.3 0.05 0.4 0.4]); hold on
pcolor(lonp,latp,los_ramp);
scatter(long,latg,sz,'k');
shading flat
colormap jet
colorbar
title('Ramp');
axis([lonmin lonmax latmin latmax]);
set(gca,'Fontsize',20);

set(gcf, 'PaperPositionMode', 'auto');

%% detrended grd files are in cm
grdwrite2(lonp,latp,los_clean,[this_track,'/ASC64_detrend.grd']);
