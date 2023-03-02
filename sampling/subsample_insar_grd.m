function subsample_insar_grd(filepath,input_grid,output_grid,grid_size)
% input_grid and output_grid should be in the same file path
% insar grid should be in the unit of cm
% usually we apply the grid size of output file is round 300 meter.
format long

    Re = 6378100;   % mean Earth radius 
    d2r = pi / 180;
    
    [lon,lat,~] = grdread2([filepath,'/',input_grid]);
    xinc = abs(lon(2) - lon(1));
    yinc = abs(lat(2) - lat(1));
    lat0 = mean(lat);
    
    xm = Re * cosd(lat0) * xinc * d2r;
    ym = Re * yinc * d2r;
    inc_m = (xm + ym)/2;
    
    Nsamp = round(grid_size / inc_m);
    new_x = Nsamp * xinc;
    new_y = Nsamp * yinc;
    
%     lon_new = linspace(lon(1),lon(end),len(lon)/Nsamp);
%     lat_new = linspace(lat(1),lat(end),len(lat)/Nsamp);
    interval = [' -I',num2str(new_x),'/',num2str(new_y),' '];
%     disp(interval);
    
    system(['gmt grdsample ',filepath,'/',input_grid,interval,' -G',filepath,'/',output_grid]);
end