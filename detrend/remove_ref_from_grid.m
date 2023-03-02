function remove_ref_from_grid(grdin,grdout,lon,lat,ref_lon,threshold)
% This function is used to subtract the value of reference point from each .grd file
% threshold in km units.

    addpath('/Users/zej011/work/Kang_tutorial/candis');
    addpath(genpath('/Users/zej011/work/zeyu/matlab'));
    
    [glon,glat,los] = grdread2(grdin);
    [mlon,mlat] = meshgrid(glon,glat);
    
    mlon = mlon(:);
    mlat = mlat(:);
    los_flat = los(:);
%     [mx,my] = utm2ll(mlon,mlat,0,1);
%     [x0,y0] = utm2ll(lon,lat,0,1);
    [mx,my] = ll2xy(mlon,mlat,ref_lon);
    [x0,y0] = ll2xy(lon,lat,ref_lon);
    
    % relative distance to the selected point (lon/lat)
    rela_x = (mx-x0)./1000;
    rela_y = (my-y0)./1000;
    radius = sqrt(rela_x.^2+rela_y.^2);
    indx = radius <= threshold;  % threshold in km

    point_value = nanmean(los_flat(indx));
    disp(point_value);
%     point_var = nanstd(los_flat(indx));
    
    los_detrend = los - point_value;
    grdwrite2(glon,glat,los_detrend,grdout);
end