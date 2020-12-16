function sign_mask_offset(this_track,offset_grid,varargin)
% apply sign mask to the noisy offset grid
% each directory has clean_left.txt and clean_right.txt inside

    [lon,lat,offset] = grdread2([this_track,'/',offset_grid]);
    [mlon,mlat] = meshgrid(lon,lat);
    clean_left = load([this_track,'/clean_left.txt']);
    clean_right = load([this_track,'/clean_right.txt']);
    tmp = mask_phase(mlon,mlat,offset,clean_left);
    offset_clean = mask_phase(mlon,mlat,tmp,clean_right);
    
    if ~isempty(varargin)
        out_grid = lower(varargin);
    else
        out_grid = 'los_clean_detrend.grd';
    end    
    grdwrite2(lon,lat,offset_clean,[this_track,'/',out_grid]);
    
    lonmin = min(lon);  lonmax = max(lon);
    latmin = min(lat);  latmax = max(lat);
    
    figure;
    subplot(1,2,1);  hold on
    pcolor(lon,lat,offset);
    shading flat
    colormap jet
    colorbar 
    axis([lonmin lonmax latmin latmax]);
    caxis([-200 200]);
    title('Raw Detrended MAI');
    set(gca,'Fontsize',20);

    subplot(1,2,2); hold on
    pcolor(lon,lat,offset_clean);
    shading flat
    colormap jet
    colorbar
    axis([lonmin lonmax latmin latmax]);
    caxis([-200 200]);
    title('Sign Masked MAI');
    set(gca,'Fontsize',20);
    
    set(gcf,'PaperPositionMode','auto');
end