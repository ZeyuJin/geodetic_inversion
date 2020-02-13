function plot_insar_detrend(this_track,insar_file,dem_file,ramp_coef)
% make sure that insar and dem file are in the same directory of this_track
    
    [lon,lat,dem] = grdread2([this_track,'/',dem_file]);
    [lon,lat,los] = grdread2([this_track,'/',insar_file]);
    [mlon,mlat] = meshgrid(lon,lat);   
    
    lon_eq = -117.5;
    lat_eq = 35.5;
    [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);
    [xsar,ysar] = utm2ll(mlon(:),mlat(:),0,1);
    xsar = (xsar - xo) ./ 1000;
    ysar = (ysar - yo) ./ 1000;
    dem = dem ./ 1000;
    
    xsar = reshape(xsar,size(mlon));
    ysar = reshape(ysar,size(mlat));
    
    if length(ramp_coef) == 4
        ramp = ramp_coef(1).*xsar + ramp_coef(2).*ysar + ramp_coef(3).*dem + ramp_coef(4);
    elseif length(ramp_coef) == 7
        ramp = ramp_coef(1).*xsar.^2 + ramp_coef(2).*ysar.^2 + ramp_coef(3).*xsar.*ysar + ...
            ramp_coef(4).*xsar + ramp_coef(5).*ysar + ramp_coef(6).*dem + ramp_coef(7);
    elseif length(ramp_coef) == 5
        ramp = ramp_coef(1).*xsar.*ysar + ramp_coef(2).*xsar + ramp_coef(3).*ysar + ramp_coef(4).*dem + ramp_coef(5);
    else
        ramp = zeros(size(los));
    end
    los_detrend = los - ramp;
    
%     h0=figure('units','normalized','outerposition',[0 0 1 1]);
%     set(h0,'renderer','painters');
    figure;
    
    subplot('Position',[0.03 0.55 0.45 0.4]);
    pcolor(xsar,ysar,los);
    shading flat
    colormap jet
    colorbar
    title('Original inteferogram (cm)');
    set(gca,'Fontsize',20);
    
    subplot('Position',[0.53 0.55 0.45 0.4]);
    pcolor(xsar,ysar,los_detrend);
    shading flat
    colormap jet
    colorbar
    title('De-trended interferogram (cm)');
    set(gca,'Fontsize',20);
    
    subplot('Position',[0.25 0.03 0.45 0.4]);
    pcolor(xsar,ysar,ramp);
    shading flat
    colormap jet
    colorbar
    title('Fitting ramp (cm)');
    set(gca,'Fontsize',20);
    
    set(gcf,'PaperPositionMode','auto');
    
    grdwrite2(lon,lat,los_detrend,[this_track,'/los_clean_detrend.grd']);
end