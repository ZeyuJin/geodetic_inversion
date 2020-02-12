function mask_insar_phase(filepath,insar_file,mask_file,wavelength,varargin)
% insar_file and mask_file should be under the same directory of filepath
% clean the InSAR data and subsample them into 300 meters resolution
% LOS : original phase (rad)
% AZO / RNG : in CM scale

    this_track = filepath;
    [lon,lat,unw] = grdread2([this_track,'/',insar_file]);
    unw_old = unw;
    [mlon,mlat] = meshgrid(lon,lat);

    % all mask?.txt saved in a common file
    fid = fopen([this_track,'/',mask_file]);
    C = textscan(fid,'%s\n');
    mask_all = C{1};
    n_mask = length(mask_all);
    fclose(fid);

    for ii = 1:n_mask
        msk_path = [this_track,'/',mask_all{ii}];
        area = load(msk_path);
        tmp = mask_phase(mlon,mlat,unw,area,'all');
        unw = tmp;
        clear tmp
    end
    unw_clean = unw;
    unw = unw_old;
    
    % convert the phase/offsets to LOS displacements (cm)
    % Phase: scale = wavelength; offsets: scale = -4pi;
    los_clean = -unw_clean*wavelength/4/pi;    

    lonmin = min(lon); lonmax = max(lon);
    latmin = min(lat); latmax = max(lat);

    figure;
    subplot('Position',[0.03 0.3 0.45 0.5]); hold on
    pcolor(lon,lat,unw);
    shading flat
    colormap jet
    colorbar
    axis([lonmin lonmax latmin latmax]);
    set(gca,'Fontsize',15);

    subplot('Position',[0.53 0.3 0.45 0.5]); hold on
    pcolor(lon,lat,unw_clean);
    shading flat
    colormap jet
    colorbar
    axis([lonmin lonmax latmin latmax]);
    set(gca,'Fontsize',15);    
    set(gcf,'PaperPositionMode','auto');
    
    out_grid = 'unwrap_clean.grd';
    grid_size = 100;    % sub-grid the file in 100 meters resolution (default)
    if ~isempty(varargin)
        for CC = 1:floor(length(varargin)/2)
            try
                switch lower(varargin{CC*2-1})
                    case 'clean_grid_name'
                        out_grid = varargin{CC*2};
                    case 'sample_grid_size'
                        grid_size = varargin{CC*2};
                end
            catch
                error('Unrecognized Keyword');
            end
        end
    end          
    grdwrite2(lon,lat,los_clean,[this_track,'/',out_grid]); 
    
    % sample the grid into 100~200 meters as the original data
    sampled_grid = 'unwrap_clean_sample.grd';
    subsample_insar_grd(this_track,out_grid,sampled_grid,grid_size);
end
