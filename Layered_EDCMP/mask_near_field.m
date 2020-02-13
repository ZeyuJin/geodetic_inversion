function mask_near_field(slip_model_mat,data_mat,fault_file,width_zone,obs_name)
% mask out extreme near-field data (within 150m in layered model)
% obs_name like 'ASC64_los' ,'DES71_rng', etc...

    d1 = load(slip_model_mat);
    slip_model_all = d1.slip_model;
    
    d2 = load(data_mat);
    los = d2.sampled_insar_data;
    xx = los(:,1);   yy = los(:,2);  zout = los(:,3);
    ve = los(:,4);   vn = los(:,5);  vz = los(:,6);    
    xx_tmp = xx;  yy_tmp = yy;  zout_tmp = zout;    % as a copy

    lon_eq = -117.5;
    lat_eq = 35.5;    
    [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);
    % read fault data
    fault_trace = load(fault_file);
    lonf = [fault_trace(:,1);fault_trace(:,3)];  
    latf = [fault_trace(:,2);fault_trace(:,4)];
    LS = length(lonf) / 2;    
    
    sz = 20;
    figure;
    subplot(1,2,1); hold on
    scatter(xx,yy,sz,zout,'filled');
    colormap jet
    colorbar
    title('Unmasked sampled data');
    set(gca,'Fontsize',20);
    for ii = 1:LS
       slon = [lonf(ii) lonf(ii+LS)];
       slat = [latf(ii) latf(ii+LS)];
       [xf,yf] = utm2ll(slon,slat,0,1);
       xf = xf - xo;
       yf = yf - yo;
       line(xf,yf,'color','black','linewidth',1.5);
    end
    
    zstart = slip_model_all(:,6);
    shallow_id = find(zstart < 0.01 & zstart > -0.01); 
    slip_model = slip_model_all(shallow_id,:);   % compute distance only using shallow patches
    NP = size(slip_model,1);
    indx = [];
    
    for ii = 1:NP
        xstart = slip_model(ii,4);
        ystart = slip_model(ii,5);
        lp = slip_model(ii,7);
        strike = slip_model(ii,9);
        
        xmid = xstart + lp/2 * sind(strike);
        ymid = ystart + lp/2 * cosd(strike);
        
        xend = xstart + lp * sind(strike);
        yend = ystart + lp * cosd(strike);
        
        dist1 = sqrt((xx - xstart).^2 + (yy - ystart).^2);
        dist2 = sqrt((xx - xmid).^2 + (yy - ymid).^2);
        dist3 = sqrt((xx - xend).^2 + (yy - yend).^2);
        
        tmp1 = find(dist1 < width_zone);
        tmp2 = find(dist2 < width_zone);
        tmp3 = find(dist3 < width_zone);
        
        tmp = union(tmp1,tmp2);
        tmp = union(tmp,tmp3);
        indx = unique([indx;tmp(:)]);
    end
    xx_tmp(indx) = [];
    yy_tmp(indx) = [];
    zout_tmp(indx) = [];
    ve(indx) = [];
    vn(indx) = [];
    vz(indx) = [];
    disp(['Dataset: ',obs_name]);
    disp(['Number of original points: ',num2str(length(xx))]);
    disp(['Number of points within ',num2str(width_zone),'m: ',num2str(length(indx))]);
    disp(['Number of points left: ',num2str(length(xx_tmp))]);
        
    sampled_insar_data = [xx_tmp,yy_tmp,zout_tmp,ve,vn,vz];
    [filepath,filename,~] = fileparts(data_mat);
    save([filepath,'/',filename,'_mask.mat'],'sampled_insar_data');

    subplot(1,2,2); hold on
    scatter(xx(indx),yy(indx),sz,zout(indx),'filled');
    colormap jet
    colorbar
    title('Masked sampled data');
    set(gca,'Fontsize',20);
    for ii = 1:LS
       slon = [lonf(ii) lonf(ii+LS)];
       slat = [latf(ii) latf(ii+LS)];
       [xf,yf] = utm2ll(slon,slat,0,1);
       xf = xf - xo;
       yf = yf - yo;
       line(xf,yf,'color','black','linewidth',1.5);
    end
    
%     % save the observation points file as the input of EDCMP
%     temp = [yy_tmp(:),xx_tmp(:)]';    % x is North in EDCMP
%     fid = fopen([filepath,'/',obs_name,'.rec'],'w');
%     fwrite(fid,temp,'real*4');
%     fclose(fid);
    
end