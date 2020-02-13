function detrend_from_residual(this_track,varargin)
    % default values
    lon_eq = -117.5;
    lat_eq = 35.5;
    fault_file = '';
    Nlook = 3;
    res_max = 20;
    ramp_type = 'bi_ramp';
    data_type = 'insar';
    model_type = [];
    
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
    
    % compute the forward model
    [xm1,ym1] = meshgrid(lon1,lat1);
    [xutm,yutm] = utm2ll(xm1(:),ym1(:),0,1);
    [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);
    xin = xutm - xo;
    yin = yutm - yo;
    xin = reshape(xin,size(xm1));
    yin = reshape(yin,size(ym1));       
    slip_model = load('fault_M7.slip');
    if strcmp(data_type,'insar')
        los_model = slip2insar_okada(xin,yin,losl,zel,znl,zul,slip_model);
    else
        los_model = slip2AZO_okada(xin,yin,losl,zel,znl,zul,slip_model);
    end
    res = losl - los_model;     % full resolution residual
   
    % mask out near field data
    mask_polygon = load('near_field_mask.txt');
%     mask_polygon = load('CSK_mask.txt');
    lonp = mask_polygon(:,1);
    latp = mask_polygon(:,2);
    [xp,yp] = utm2ll(lonp,latp,0,1);
    xv = xp - xo;
    yv = yp - yo;
    out = ~inpolygon(xin,yin,xv,yv);
    xout = xin(out) ./ 1000;
    yout = yin(out) ./ 1000;
    res_out = res(out);  
    
    % convert to km scale for plot and de-trend
    xin_pl = xin ./ 1000;
    yin_pl = yin ./ 1000;
    xv_pl = xv ./ 1000;
    yv_pl = yv ./ 1000;
   
    % detrend from the residual
    if strcmp(ramp_type,'bi_ramp')
       pfit = fit_ramp(res_out,xout,yout);
       ramp = pfit(1).*xin_pl + pfit(2).*yin_pl + pfit(3);
    elseif strcmp(ramp_type,'qu_ramp_5')
       pfit = fit_ramp_hyper(res_out,xout,yout);
       ramp = pfit(1).*xin_pl.*yin_pl + pfit(2).*xin_pl + pfit(3).*yin_pl + pfit(4);
    elseif strcmp(ramp_type,'qu_ramp_7')
       pfit = fit_ramp_quadratic(res_out,xout,yout);
       ramp = pfit(1).*xin_pl.^2 + pfit(2).*yin_pl.^2 + pfit(3).*xin_pl.*yin_pl + pfit(4).*xin_pl + pfit(5).*yin_pl + pfit(6);
    else
       ramp = zeros(size(res));
       pfit = 0;
    end  
    res_second = res - ramp;
    
    figure;
    subplot('Position',[0.03 0.55 0.45 0.4]); hold on
    pcolor(xin/1000,yin/1000,res);
    shading flat
    colormap jet
    colorbar
    title('Residual from misfit (cm)');
    set(gca,'Fontsize',20);
    caxis([-res_max res_max]);
%     caxis([-150 150]);
    if ~isempty(fault_file)
       for ii = 1:LS
          slon = [lonf(ii) lonf(ii+LS)];
          slat = [latf(ii) latf(ii+LS)];
          [xx,yy] = utm2ll(slon,slat,0,1);
          xs = (xx - xo) ./ 1000;
          ys = (yy - yo) ./ 1000;
          line(xs,ys,'color','black','linewidth',1.5);
       end
    end
%     plot(xv_pl,yv_pl,'k','linewidth',1);  % plot the mask area
    
    subplot('Position',[0.53 0.55 0.45 0.4]); hold on
    pcolor(xin/1000,yin/1000,ramp);
    shading flat
    colormap jet
    colorbar
    title('Ramp from residual (cm)');
    set(gca,'Fontsize',20);
%     caxis([-res_max res_max]);
%     caxis([-150 150]);
    if ~isempty(fault_file)
       for ii = 1:LS
          slon = [lonf(ii) lonf(ii+LS)];
          slat = [latf(ii) latf(ii+LS)];
          [xx,yy] = utm2ll(slon,slat,0,1);
          xs = (xx - xo) ./ 1000;
          ys = (yy - yo) ./ 1000;
          line(xs,ys,'color','black','linewidth',1.5);
       end
    end
%     plot(xv_pl,yv_pl,'k','linewidth',1);
    
    subplot('Position',[0.25 0.03 0.45 0.4]); hold on
    pcolor(xin/1000,yin/1000,res_second);
    shading flat
    colormap jet
    colorbar
    title('Final residual (cm)');
    set(gca,'Fontsize',20);
    caxis([-res_max res_max]);
    if ~isempty(fault_file)
       for ii = 1:LS
          slon = [lonf(ii) lonf(ii+LS)];
          slat = [latf(ii) latf(ii+LS)];
          [xx,yy] = utm2ll(slon,slat,0,1);
          xs = (xx - xo) ./ 1000;
          ys = (yy - yo) ./ 1000;
          line(xs,ys,'color','black','linewidth',1.5);
       end
    end
%     plot(xv_pl,yv_pl,'k','linewidth',1);    
    set(gcf,'PaperPositionMode','auto');
    
    % save the final residual and ramp coef
    save([this_track,'/residual_fit.mat'],'pfit');
    grdwrite2(lon1,lat1,res_second,[this_track,'/los_residual_',model_type,'.grd']);
end