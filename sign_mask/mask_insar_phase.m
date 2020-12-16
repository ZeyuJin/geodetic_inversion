function mask_insar_phase(filepath,insar_file,mask_file,wavelength,varargin)
% insar_file and mask_file should be under the same directory of filepath
% clean the InSAR data and subsample them into 300 meters resolution
% LOS : original phase (rad)
% AZO / RNG : in CM scale

    mask_type = 'allout';
    out_grid = 'unwrap_clean.grd';
    grid_size = 100;    % sub-grid the file in 100 meters resolution (default)   
    los_max = 6;  % in cm
    detrend = 0;
    
    if ~isempty(varargin)
        for CC = 1:floor(length(varargin)/2)
            try
                switch lower(varargin{CC*2-1})
                    case 'clean_grid_name'
                        out_grid = varargin{CC*2};
                    case 'sample_grid_size'
                        grid_size = varargin{CC*2};
                    case 'mask_type'
                        mask_type = varargin{CC*2};
                    case 'los_max'
                        los_max = varargin{CC*2};
                    case 'detrend'
                        detrend = varargin{CC*2};
                end
            catch
                error('Unrecognized Keyword');
            end
        end
    end   

    this_track = filepath;
    [lon,lat,unw] = grdread2([this_track,'/',insar_file]);
    unw_old = unw;
    [mlon,mlat] = meshgrid(lon,lat);
    
%     % if detrend the phase with topo
%     if detrend
%         if isfile([this_track,'/dem_samp.grd'])
%             [~,~,topo] = grdread2([this_track,'/dem_samp.grd']);
%             topo = topo./1000;
%             pout = fit_ramp_topo(unw,mlon,mlat,topo);
%             unw = unw - (pout(1).*mlon + pout(2).*mlat + pout(3).*topo + pout(4));
%         else
%             disp('No dem_samp.grd found! Do not detrend the phase!');
%         end          
%     end

    % all mask?.txt saved in a common file
    fid = fopen([this_track,'/',mask_file]);
    C = textscan(fid,'%s\n');
    mask_all = C{1};
    n_mask = length(mask_all);
    fclose(fid);

%     unw_clean = unw;
    for ii = 1:n_mask
        msk_path = [this_track,'/',mask_all{ii}];
        area = load(msk_path);
        tmp = mask_phase(mlon,mlat,unw,area,mask_type);
        unw = tmp;
        clear tmp
    end
    unw_clean = unw;
%     unw = unw_old;
    
    % if detrend the phase with topo
    ramp = zeros(size(unw_clean));
    if detrend
        if isfile([this_track,'/dem_samp.grd'])
            [~,~,topo] = grdread2([this_track,'/dem_samp.grd']);
            topo = topo./1000;
            pout = fit_ramp_topo(unw_clean(:),mlon(:),mlat(:),topo(:));
            ramp = pout(1).*mlon + pout(2).*mlat + pout(3).*topo + pout(4);
            unw_clean = unw_clean - ramp;
        else
            disp('No dem_samp.grd found! Do not detrend the phase!');
        end          
    end
    
    % convert the phase/offsets to LOS displacements (cm)
    % Phase: scale = wavelength; offsets: scale = -4pi;
    los_clean = -unw_clean*wavelength/4/pi;   
    grdwrite2(lon,lat,los_clean,[this_track,'/',out_grid]);

    lonmin = min(lon); lonmax = max(lon);
    latmin = min(lat); latmax = max(lat);

    figure;
%     subplot('Position',[0.03 0.3 0.45 0.5]); hold on
    subplot(2,2,1); hold on
    pcolor(lon,lat,unw_old);
    shading flat
    colormap jet
    colorbar
    axis([lonmin lonmax latmin latmax]);
    caxis([-los_max los_max]);
    title('Original LOS');
    set(gca,'Fontsize',15);

%     subplot('Position',[0.53 0.3 0.45 0.5]); hold on
    subplot(2,2,2); hold on
    pcolor(lon,lat,unw);
    shading flat
    colormap jet
%     colorbar
    axis([lonmin lonmax latmin latmax]);
    caxis([-los_max los_max]);
    title('Masked LOS');
    set(gca,'Fontsize',15);    
    
    subplot(2,2,3); hold on
    pcolor(lon,lat,unw_clean);
    shading flat
    colormap jet
    colorbar
    axis([lonmin lonmax latmin latmax]);
    caxis([-los_max los_max]);
    title('Masked and Detrended LOS');
    set(gca,'Fontsize',15); 
    
    subplot(2,2,4); hold on
    pcolor(lon,lat,ramp);
    shading flat
    colormap jet
%     colorbar
    axis([lonmin lonmax latmin latmax]);
    caxis([-los_max los_max]);
    title('Detrended Ramp');
    set(gca,'Fontsize',15);
        
    set(gcf,'PaperPositionMode','auto');     
    
    % sample the grid into 100~200 meters as the original data
    sampled_grid = 'unwrap_clean_sample.grd';
    subsample_insar_grd(this_track,out_grid,sampled_grid,grid_size);
end
