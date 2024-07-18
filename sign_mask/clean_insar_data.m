wavelength_C = 0.0554658*100;  % wavelength of C-band
wavelength_L = 0.242452*100;   % wavelength of L-band
scale = -4*pi;                 % for offsets data

% this_track = '/Users/zej011/coseismic/ASC100/LOS';
% insar_file = 'unwrap_ll.grd';
% mask_file = 'mask_txt';
% mask_insar_phase(this_track,insar_file,mask_file,wavelength_C);

% this_track = '/Users/zej011/coseismic/ASC100/LOS2';
% insar_file = 'unwrap_ll.grd';
% mask_file = 'mask_txt';
% mask_insar_phase(this_track,insar_file,mask_file,wavelength_C,'los_max',80);

% this_track = '/Users/zej011/coseismic/DES5/LOS';
% insar_file = 'unwrap_ll.grd';
% mask_file = 'mask_txt';
% mask_insar_phase(this_track,insar_file,mask_file,wavelength_C);

% this_track = '/Users/zej011/coseismic/DES5/LOS2';
% insar_file = 'unwrap_ll.grd';
% mask_file = 'mask_txt';
% mask_insar_phase(this_track,insar_file,mask_file,wavelength_C,'los_max',80);

this_track = '/Users/zej011/coseismic/DES5/LOS3';
insar_file = 'unwrap_ll.grd';
mask_file = 'mask_txt';
mask_insar_phase(this_track,insar_file,mask_file,wavelength_C,'los_max',80,'detrend',1,'mask');

% this_track = '/Users/zej011/coseismic/ASC100/offsets';
% insar_file = 'rng_off_filt_ll.grd';
% mask_file = 'mask_txt';
% mask_insar_phase(this_track,insar_file,mask_file,scale);

% this_track = '/Users/zej011/coseismic/DES5/offsets';
% insar_file = 'rng_off_filt_ll.grd';
% mask_file = 'mask_txt';
% mask_insar_phase(this_track,insar_file,mask_file,scale);

% this_track = '/Users/zej011/coseismic/ALOS2_SCAN';
% insar_file = 'unwrap_corrected_ll.grd';
% mask_file = 'mask_txt';
% mask_insar_phase(this_track,insar_file,mask_file,wavelength_L);

% this_track = '/Users/zej011/coseismic/ALOS2_SCAN2';
% insar_file = 'merge_unwrap_ll.grd';
% mask_file = 'mask_txt';
% mask_insar_phase(this_track,insar_file,mask_file,wavelength_L,'los_max',60);

% this_track = '/Users/zej011/coseismic/ALOS2_stripe/range';
% insar_file = 'rng_offset_ll.grd';
% mask_file = 'mask_txt';
% mask_insar_phase(this_track,insar_file,mask_file,scale);

% this_track = '/Users/zej011/coseismic/ALOS2_stripe/azimuth';
% insar_file = 'azi_offset_ll.grd';
% mask_file = 'mask_txt';
% mask_insar_phase(this_track,insar_file,mask_file,scale);

% this_track = '/Users/zej011/coseismic/ALOS2_stripe/LOS';
% insar_file = 'unwrap_ll.grd';
% mask_file = 'mask_txt';
% mask_insar_phase(this_track,insar_file,mask_file,wavelength_L);

% [lon,lat,raw] = grdread2([this_track,'/',insar_file]);
% [mlon,mlat] = meshgrid(lon,lat);
% clean_left = load([this_track,'/clean_left.txt']);
% clean_right = load([this_track,'/clean_right.txt']);
% raw = raw - 20;    % shift the mean value
% tmp = mask_phase(mlon,mlat,raw,clean_left);
% raw_clean = mask_phase(mlon,mlat,tmp,clean_right);
% raw_clean = raw_clean + 20;
% raw = raw + 20;
% clear tmp
% 
% grdwrite2(lon,lat,raw_clean,[this_track,'/ASC_rng_ll_clean.grd']);
% insar_file = 'ASC_rng_ll_clean.grd';
% scale = -4*pi;       % because it is the offset not phase
% mask_insar_phase(this_track,insar_file,mask_file,scale);
% 
% this_track = '/Users/zej011/Ridgecrest/data_resample/DES71/offsets';
% insar_file = 'DES_rng_ll_finer.grd';
% mask_file = 'maskfile';
% mask_insar_phase(this_track,insar_file,mask_file,scale);
% 
% [lon,lat,raw] = grdread2([this_track,'/',insar_file]);
% [mlon,mlat] = meshgrid(lon,lat);
% clean_left = load([this_track,'/clean_left.txt']);
% clean_right = load([this_track,'/clean_right.txt']);
% tmp = mask_phase(mlon,mlat,raw,clean_left);
% raw_clean = mask_phase(mlon,mlat,tmp,clean_right);
% clear tmp
% 
% grdwrite2(lon,lat,raw_clean,[this_track,'/DES_rng_ll_clean.grd']);
% insar_file = 'DES_rng_ll_clean.grd';
% scale = -4*pi;       % because it is the offset not phase
% mask_insar_phase(this_track,insar_file,mask_file,scale); 

% this_track = '/Users/zej011/Ridgecrest/data_resample/Cosmo_Skymed/ASC_offsets';
% insar_file = 'ASC_azi_filt.grd';
% mask_file = 'maskfile';
% mask_insar_phase(this_track,insar_file,mask_file,scale);

% [lon,lat,ASC_azi] = grdread2([this_track,'/',insar_file]);
% [mlon,mlat] = meshgrid(lon,lat);
% area = load([this_track,'/ASC_mask.txt']);
% out = ~inpolygon(mlon,mlat,area(:,1),area(:,2));
% ASC_azi = ASC_azi * 100;
% ASC_azi_mask = ASC_azi;
% ASC_azi_mask(out) = NaN;
% 
% mask_left = load([this_track,'/ASC_mask_left']);
% mask_right = load([this_track,'/ASC_mask_right']);
% tmp = mask_phase(mlon,mlat,ASC_azi_mask,mask_left);
% ASC_azi_clean_mask = mask_phase(mlon,mlat,tmp,mask_right);
% clear tmp
% grdwrite2(lon,lat,ASC_azi_clean_mask,[this_track,'/unwrap_clean_sample.grd']);
% 
% lonmin = min(lon); lonmax = max(lon);
% latmin = min(lat); latmax = max(lat);
% 
% figure;
% subplot('Position',[0.035 0.55 0.45 0.4]);  hold on
% pcolor(lon,lat,ASC_azi);
% shading flat
% colormap jet
% colorbar 
% axis([lonmin lonmax latmin latmax]);
% title('Azimuth offsets (ASC)');
% set(gca,'Fontsize',20);
% caxis([-200 200]);
% 
% subplot('Position',[0.535 0.55 0.45 0.4]); hold on
% pcolor(lon,lat,ASC_azi_mask);
% shading flat
% colormap jet
% colorbar
% axis([lonmin lonmax latmin latmax]);
% title('Masked Azimuth offsets (ASC)');
% set(gca,'Fontsize',20);
% caxis([-200 200]);
% 
% subplot('Position',[0.25 0.03 0.45 0.4]); hold on
% pcolor(lon,lat,ASC_azi_clean_mask);
% shading flat
% colormap jet
% colorbar
% axis([lonmin lonmax latmin latmax]);
% title('Cleaned Masked Azimuth offsets (ASC)');
% set(gca,'Fontsize',20);
% caxis([-200 200]);
% set(gcf,'PaperPositionMode','auto');

% this_track = '/Users/zej011/Ridgecrest/data_resample/Cosmo_Skymed/DES_offsets';
% insar_file = 'DES_azi_filt.grd';
% mask_file = 'maskfile';
% mask_insar_phase(this_track,insar_file,mask_file,scale);

% [lon,lat,DES_azi] = grdread2([this_track,'/',insar_file]);
% [mlon,mlat] = meshgrid(lon,lat);
% area = load([this_track,'/DES_mask.txt']);
% out = ~inpolygon(mlon,mlat,area(:,1),area(:,2));
% DES_azi = DES_azi * 100;
% DES_azi_mask = DES_azi;
% DES_azi_mask(out) = NaN;
% 
% mask_left = load([this_track,'/DES_mask_left']);
% mask_right = load([this_track,'/DES_mask_right']);
% tmp = mask_phase(mlon,mlat,DES_azi_mask,mask_left);
% DES_azi_clean_mask = mask_phase(mlon,mlat,tmp,mask_right);
% clear tmp
% grdwrite2(lon,lat,DES_azi_clean_mask,[this_track,'/unwrap_clean_sample.grd']);
% 
% lonmin = min(lon); lonmax = max(lon);
% latmin = min(lat); latmax = max(lat);
% 
% figure;
% subplot('Position',[0.035 0.55 0.45 0.4]);  hold on
% pcolor(lon,lat,DES_azi);
% shading flat
% colormap jet
% colorbar 
% axis([lonmin lonmax latmin latmax]);
% title('Azimuth offsets (DES)');
% set(gca,'Fontsize',20);
% caxis([-200 200]);
% 
% subplot('Position',[0.535 0.55 0.45 0.4]); hold on
% pcolor(lon,lat,DES_azi_mask);
% shading flat
% colormap jet
% colorbar
% axis([lonmin lonmax latmin latmax]);
% title('Masked Azimuth offsets (DES)');
% set(gca,'Fontsize',20);
% caxis([-200 200]);
% 
% subplot('Position',[0.25 0.03 0.45 0.4]); hold on
% pcolor(lon,lat,DES_azi_clean_mask);
% shading flat
% colormap jet
% colorbar
% axis([lonmin lonmax latmin latmax]);
% title('Cleaned Masked Azimuth offsets (DES)');
% set(gca,'Fontsize',20);
% caxis([-200 200]);
% set(gcf,'PaperPositionMode','auto');
