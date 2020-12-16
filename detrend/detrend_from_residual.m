function detrend_from_residual(this_track,varargin)
% there is something wrong when I computed the forward model
% using layered solutions (call slip2insar_okada.m before)
% add module "slip2disp_edcmp.m"
% fixed by Zeyu Jin on Mar 25th, 2020

    % default values
    lon_eq = -117.5;
    lat_eq = 35.5;
    ref_lon = lon_eq;
    fault_file = '';
    Nlook = 3;
    res_max = 20;
    ramp_type = 'bi_ramp';
    data_type = 'insar';
    model_type = 'okada';
    near_mask = 0;
    
    if ~isempty(varargin)
       for CC = 1:floor(length(varargin)/2)
           try
               switch lower(varargin{CC*2-1})
                   case 'lonc'
                       lon_eq = varargin{CC*2};
                   case 'latc'
                       lat_eq = varargin{CC*2};
                   case 'misfit_range'
                       res_max = varargin{CC*2};
                   case 'ramp_type'
                       ramp_type = varargin{CC*2};
                   case 'dec'
                       Nlook = varargin{CC*2};
                   case 'fault'
                       fault_file = varargin{CC*2};
                   case 'data_type'
                       data_type = varargin{CC*2};
                   case 'model_type'
                       model_type = varargin{CC*2};
                   case 'ref_lon'
                       ref_lon = varargin{CC*2};
                   case 'near_mask'
                       near_mask = varargin{CC*2};
               end
           catch
               error('Unrecognized Keyword');
           end
       end
    end
   
   % read fault data
   if ~isempty(fault_file)
       fault_trace = load(fault_file);
       lonf = [fault_trace(:,1);fault_trace(:,3)];  
       latf = [fault_trace(:,2);fault_trace(:,4)];
       LS = length(lonf) / 2;
   end
   
    % read full resolution data
    [x1,y1,losin]=grdread2([this_track,'/','los_ll_low','.grd']);   % in the unit of cm
    [~,~,ze]=grdread2([this_track,'/','look_e_low','.grd']);
    [~,~,zn]=grdread2([this_track,'/','look_n_low','.grd']);
    [~,~,zu]=grdread2([this_track,'/','look_u_low','.grd']);
    
    % multi-look to reduce the computation time
    [lon1,lat1,losl] = multi_look(x1,y1,losin,Nlook,Nlook);
    [~,~,zel] = multi_look(x1,y1,ze,Nlook,Nlook);
    [~,~,znl] = multi_look(x1,y1,zn,Nlook,Nlook);
    [~,~,zul] = multi_look(x1,y1,zu,Nlook,Nlook);
    
    grdwrite2(lon1,lat1,zel,[this_track,'/ze_low.grd']);
    grdwrite2(lon1,lat1,znl,[this_track,'/zn_low.grd']);
    grdwrite2(lon1,lat1,zul,[this_track,'/zu_low.grd']);
    
    % compute the forward model
    [xm1,ym1] = meshgrid(lon1,lat1);
%     [xutm,yutm] = utm2ll(xm1(:),ym1(:),0,1);
%     [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);
    [xutm,yutm] = ll2xy(xm1(:),ym1(:),ref_lon);
    [xo,yo] = ll2xy(lon_eq,lat_eq,ref_lon);
    xin = xutm - xo;
    yin = yutm - yo;
    xin = reshape(xin,size(xm1));
    yin = reshape(yin,size(ym1));       
%     slip_model = load('fault_M7.slip');
%     slip_model = load('fault_M5.slip');
%     tmp = load('resample/4_segments/homo_4data.mat');
%     tmp = load('resample/misfit_include_MAI/homo_better_data.mat');
    tmp = load('resample/two_segments/afterslip_strike.mat');
    slip_model = tmp.slip_model;
    
    
    if strcmp(model_type,'okada')
        if strcmp(data_type,'insar')
            los_model = slip2insar_okada(xin,yin,losl,zel,znl,zul,slip_model);
        else
            los_model = slip2AZO_okada(xin,yin,losl,zel,znl,zul,slip_model);
        end
    else
        los_model = slip2disp_edcmp(this_track,xin,yin,losl,zel,znl,zul,slip_model,data_type);
    end
    res = losl - los_model;     % full resolution residual
   
    if near_mask
        % mask out near field data
        mask_polygon = load([this_track,'/near_field_mask.txt']);
        lonp = mask_polygon(:,1);
        latp = mask_polygon(:,2);
        [xp,yp] = ll2xy(lonp,latp,ref_lon);
        xv = xp - xo;
        yv = yp - yo;
        out = ~inpolygon(xin,yin,xv,yv);
        xout = xin(out) ./ 1000;
        yout = yin(out) ./ 1000;
        res_out = res(out); 
    else
        xout = xin;
        yout = yin;
        res_out = res;
    end
    
    % convert to km scale for plot and de-trend
    xin_pl = xin ./ 1000;
    yin_pl = yin ./ 1000;
    
    % detrend resampled data
%     if strcmp(model_type,'okada')
%         data = load([this_track,'/los_samp1.mat']);
%     else
%         data = load([this_track,'/los_samp1_mask.mat']);
%     end
    data = load([this_track,'/los_samp1.mat']);
    los_samp = data.sampled_insar_data(:,3);
    XS = data.sampled_insar_data(:,1) / 1000;
    YS = data.sampled_insar_data(:,2) / 1000;
    look_angle = data.sampled_insar_data(:,4:6);    
   
    % detrend from the residual
    if strcmp(ramp_type,'bi_ramp')
       pfit = fit_ramp(res_out,xout,yout);
       ramp = pfit(1).*xin_pl + pfit(2).*yin_pl + pfit(3);
       RS = pfit(1).*XS + pfit(2).*YS + pfit(3);     
    elseif strcmp(ramp_type,'qu_ramp_5')
       pfit = fit_ramp_hyper(res_out,xout,yout);
       ramp = pfit(1).*xin_pl.*yin_pl + pfit(2).*xin_pl + pfit(3).*yin_pl + pfit(4);
       RS = pfit(1).*XS.*YS + pfit(2).*XS + pfit(3).*YS + pfit(4);
    elseif strcmp(ramp_type,'qu_ramp_7')
       pfit = fit_ramp_quadratic(res_out,xout,yout);
       ramp = pfit(1).*xin_pl.^2 + pfit(2).*yin_pl.^2 + pfit(3).*xin_pl.*yin_pl + pfit(4).*xin_pl + pfit(5).*yin_pl + pfit(6);
       RS = pfit(1).*XS.^2 + pfit(2).*YS.^2 + pfit(3).*XS.*YS + pfit(4).*XS + pfit(5).*YS + pfit(6);
    else
       ramp = zeros(size(res));
       RS = zeros(size(los_samp));
       pfit = 0;
    end  
    res_second = res - ramp;
    data_detr = losl - ramp;
    los_samp_detrend = los_samp - RS;
    sampled_insar_data = double([XS*1000,YS*1000,los_samp_detrend,look_angle]);

    if strcmp(model_type,'okada')
        save_data_name = 'los_samp1_detrend.mat';
    else
        save_data_name = ['los_samp1_detrend_',model_type,'.mat'];
    end
    save([this_track,'/',save_data_name],'sampled_insar_data');
    
    figure;
    subplot('Position',[0.03 0.55 0.45 0.4]); hold on
%     pcolor(xin/1000,yin/1000,res);
    pcolor(xin/1000,yin/1000,losl);
    shading flat
    colormap jet
    colorbar
%     title('Residual from misfit (cm)');
%     plot(xv/1000,yv/1000,'linewidth',1.5,'color','m');
    title('Original along-track interferometry (cm)');
    set(gca,'Fontsize',20);
    caxis([-res_max res_max]);
%     caxis([-150 150]);
    if ~isempty(fault_file)
       for ii = 1:LS
          slon = [lonf(ii) lonf(ii+LS)];
          slat = [latf(ii) latf(ii+LS)];
%           [xx,yy] = utm2ll(slon,slat,0,1);
          [xx,yy] = ll2xy(slon,slat,ref_lon);
          xs = (xx - xo) ./ 1000;
          ys = (yy - yo) ./ 1000;
          line(xs,ys,'color','black','linewidth',1);
       end
    end
%     plot(xv_pl,yv_pl,'k','linewidth',1);  % plot the mask area
    
    subplot('Position',[0.53 0.55 0.45 0.4]); hold on
%     pcolor(xin/1000,yin/1000,ramp);
    pcolor(xin/1000,yin/1000,los_model);
    shading flat
    colormap jet
    colorbar
%     title('Ramp from residual (cm)');
    title('Model Prediction (without using MAI)');
    set(gca,'Fontsize',20);
    caxis([-res_max res_max]);
%     caxis([-150 150]);
    if ~isempty(fault_file)
       for ii = 1:LS
          slon = [lonf(ii) lonf(ii+LS)];
          slat = [latf(ii) latf(ii+LS)];
%           [xx,yy] = utm2ll(slon,slat,0,1);
          [xx,yy] = ll2xy(slon,slat,ref_lon);
          xs = (xx - xo) ./ 1000;
          ys = (yy - yo) ./ 1000;
          line(xs,ys,'color','black','linewidth',1);
       end
    end
%     plot(xv_pl,yv_pl,'k','linewidth',1);
    
    subplot('Position',[0.03 0.03 0.45 0.4]); hold on
%     pcolor(xin/1000,yin/1000,res_second);
    pcolor(xin/1000,yin/1000,res);
    shading flat
    colormap jet
    colorbar
%     title('Final residual (cm)');
    title('Residual (data - model)');
    set(gca,'Fontsize',20);
    caxis([-res_max res_max]);
    if ~isempty(fault_file)
       for ii = 1:LS
          slon = [lonf(ii) lonf(ii+LS)];
          slat = [latf(ii) latf(ii+LS)];
%           [xx,yy] = utm2ll(slon,slat,0,1);
          [xx,yy] = ll2xy(slon,slat,ref_lon);
          xs = (xx - xo) ./ 1000;
          ys = (yy - yo) ./ 1000;
          line(xs,ys,'color','black','linewidth',1);
       end
    end
    
    subplot('Position',[0.53 0.03 0.45 0.4]); hold on
%     pcolor(xin/1000,yin/1000,res_second);
%     pcolor(xin/1000,yin/1000,data_detr);
    pcolor(xin/1000, yin/1000, ramp);
    shading flat
    colormap jet
    colorbar
%     title('Final residual (cm)');
    title('Ramp');
    set(gca,'Fontsize',20);
%     caxis([-res_max res_max]);
    if ~isempty(fault_file)
       for ii = 1:LS
          slon = [lonf(ii) lonf(ii+LS)];
          slat = [latf(ii) latf(ii+LS)];
%           [xx,yy] = utm2ll(slon,slat,0,1);
          [xx,yy] = ll2xy(slon,slat,ref_lon);
          xs = (xx - xo) ./ 1000;
          ys = (yy - yo) ./ 1000;
          line(xs,ys,'color','black','linewidth',1);
       end
    end
    
%     plot(xv_pl,yv_pl,'k','linewidth',1);    
    set(gcf,'PaperPositionMode','auto');
       
    % save the final residual and ramp coef
%     save([this_track,'/residual_fit.mat'],'pfit');
    grdwrite2(lon1,lat1,res_second,[this_track,'/los_residual_',model_type,'.grd']);
    grdwrite2(lon1,lat1,data_detr,[this_track,'/los_ll_low_detrend_',model_type,'.grd']);
end