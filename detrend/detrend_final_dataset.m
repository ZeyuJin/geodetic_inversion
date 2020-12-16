function detrend_final_dataset(this_track,varargin)
    % default values
    lon_eq = -117.5;
    lat_eq = 35.5;
    ref_lon = lon_eq;
    fault_file = '';
    res_max = 20;
    axis_range = [-100 100 -100 100];
    threshold = 50;
    model_type = 'okada';
    
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
                   case 'fault'
                       fault_file = varargin{CC*2};
                   case 'axis_range'
                       axis_range = varargin{CC*2};
                       if length(axis_range) ~= 4
                           error('Something wrong with axis range');
                       end
                   case 'threshold'
                       threshold = varargin{CC*2};
                   case 'model_type'
                       model_type = varargin{CC*2};
                   case 'ref_lon'
                       ref_lon = varargin{CC*2};
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
   
   % load full resoluton residual and pfit
   tmp = load([this_track,'/residual_fit.mat']);
   pfit = tmp.pfit;
   
   % read full resolution data and detrend it
   [x1,y1,losin]=grdread2([this_track,'/los_ll_low.grd']);   % in the unit of cm
   [xm1,ym1] = meshgrid(x1,y1);
%    [xutm,yutm] = utm2ll(xm1(:),ym1(:),0,1);
%    [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);
   [xutm,yutm] = ll2xy(xm1(:),ym1(:),ref_lon);
   [xo,yo] = ll2xy(lon_eq,lat_eq,ref_lon);
    xin = (xutm - xo) ./ 1000;
    yin = (yutm - yo) ./ 1000;
    xin = reshape(xin,size(xm1));
    yin = reshape(yin,size(ym1));
    
    if length(pfit) == 3 && norm(pfit,1) ~= 0
        ramp = pfit(1).*xin + pfit(2).*yin + pfit(3);
    elseif length(pfit) == 4
        ramp = pfit(1).*xin.*yin + pfit(2).*xin + pfit(3).*yin + pfit(4);
    elseif length(pfit) == 6
        ramp = pfit(1).*xin.^2 + pfit(2).*yin.^2 + pfit(3).*xin.*yin + pfit(4).*xin + pfit(5).*yin + pfit(6);
    else
        ramp = zeros(size(xin));
    end
    los_detrend = losin - ramp;
    grdwrite2(x1,y1,los_detrend,[this_track,'/los_ll_low_detrend_',model_type,'.grd']);   % write into grid
    
%     % plot for test
%     subplot('Position',[0.03 0.55 0.45 0.4]); hold on
%     pcolor(xin,yin,losin);
%     shading flat
%     colormap jet
%     colorbar
%     title('Residual from misfit (cm)');
%     set(gca,'Fontsize',20);    
%     
%     subplot('Position',[0.53 0.55 0.45 0.4]); hold on
%     pcolor(xin,yin,los_detrend);
%     shading flat
%     colormap jet
%     colorbar
%     title('Ramp from residual (cm)');
%     set(gca,'Fontsize',20);
%     
%     subplot('Position',[0.25 0.03 0.45 0.4]); hold on
%     pcolor(xin,yin,ramp);
%     shading flat
%     colormap jet
%     colorbar
%     title('Final residual (cm)');
%     set(gca,'Fontsize',20);
%     caxis([-res_max res_max]);
    
%     % detrend resampled data
%     if strcmp(model_type,'okada')
%         data = load([this_track,'/los_samp3.mat']);
%     else
%         data = load([this_track,'/los_samp3_mask.mat']);
%     end
%     los_samp = data.sampled_insar_data(:,3);
%     XS = data.sampled_insar_data(:,1) / 1000;
%     YS = data.sampled_insar_data(:,2) / 1000;
%     look_angle = data.sampled_insar_data(:,4:6);
%     
%     if length(pfit) == 3 && norm(pfit,1) ~= 0
%         RS = pfit(1).*XS + pfit(2).*YS + pfit(3);
%     elseif length(pfit) == 4
%         RS = pfit(1).*XS.*YS + pfit(2).*XS + pfit(3).*YS + pfit(4);
%     elseif length(pfit) == 6
%         RS = pfit(1).*XS.^2 + pfit(2).*YS.^2 + pfit(3).*XS.*YS + pfit(4).*XS + pfit(5).*YS + pfit(6);
%     else
%         RS = zeros(size(los_samp));
%     end
%     los_samp_detrend = los_samp - RS;
%     sampled_insar_data = double([XS*1000,YS*1000,los_samp_detrend,look_angle]);
%     if strcmp(model_type,'okada')
%         save_data_name = 'los_samp3_detrend.mat';
%         model_name = 'los_model.mat';
%     else
%         save_data_name = ['los_samp3_detrend_',model_type,'.mat'];
%         model_name = 'los_model_layered.mat';
%     end    
%     save([this_track,'/',save_data_name],'sampled_insar_data');
%     
%     model = load([this_track,'/',model_name]);
%     los_model = model.sampled_model(:,3);
%     los_resid = los_samp_detrend - los_model;
%     
%     % delete some points which sampled on the boundary
%     indx = find(los_resid > threshold);
%     los_samp_detrend(indx) = NaN;
%     los_model(indx) = NaN;
%     
%     if strcmp(model_type,'okada')
%         resi_grd = 'los_residual_A_okada.grd';
%     else
%         resi_grd = 'los_residual_A_layer.grd';
%     end    
%     [xfull,yfull,los_resid_full] = grdread2([this_track,'/',resi_grd]);
%     [xmf,ymf] = meshgrid(xfull,yfull);
% %     [XU,YU] = utm2ll(xmf(:),ymf(:),0,1);
%     [XU,YU] = ll2xy(xmf(:),ymf(:),ref_lon);
%     XF = (XU - xo) ./ 1000;
%     YF = (YU - yo) ./ 1000;
%     XF = reshape(XF,size(xmf));
%     YF = reshape(YF,size(ymf));
%     los_resid_full(los_resid_full > threshold) = NaN;
%     
%     % plot the resampled data, model and full residual
%     sz = 20;
%     figure;
%     subplot('Position',[0.03 0.55 0.45 0.4]); hold on
%     scatter(XS,YS,sz,los_samp_detrend,'filled');
%     shading flat
%     colormap jet
%     colorbar
%     title('Sampled data (cm)');
%     set(gca,'Fontsize',20); 
%     axis(axis_range);
% %     caxis([-80 80]);
%     caxis([-2 2]);
%     if ~isempty(fault_file)
%        for ii = 1:LS
%           slon = [lonf(ii) lonf(ii+LS)];
%           slat = [latf(ii) latf(ii+LS)];
% %           [xx,yy] = utm2ll(slon,slat,0,1);
%           [xx,yy] = ll2xy(slon,slat,ref_lon);
%           xs = (xx - xo) ./ 1000;
%           ys = (yy - yo) ./ 1000;
%           line(xs,ys,'color','black','linewidth',1.5);
%        end
%     end
%     
%     subplot('Position',[0.53 0.55 0.45 0.4]); hold on
%     scatter(XS,YS,sz,los_model,'filled');
%     shading flat
%     colormap jet
%     colorbar
%     title('Sampled model (cm)');
%     set(gca,'Fontsize',20);
% %     caxis([-80 80]);
%     caxis([-2 2]);
%     axis(axis_range);
%     if ~isempty(fault_file)
%        for ii = 1:LS
%           slon = [lonf(ii) lonf(ii+LS)];
%           slat = [latf(ii) latf(ii+LS)];
% %           [xx,yy] = utm2ll(slon,slat,0,1);
%           [xx,yy] = ll2xy(slon,slat,ref_lon);
%           xs = (xx - xo) ./ 1000;
%           ys = (yy - yo) ./ 1000;
%           line(xs,ys,'color','black','linewidth',1.5);
%        end
%     end
%     
%     subplot('Position',[0.25 0.03 0.45 0.4]); hold on
%     pcolor(XF,YF,los_resid_full);
%     shading flat
%     colormap jet
%     colorbar
%     title('Full residual (cm)');
%     set(gca,'Fontsize',20);
%     caxis([-res_max res_max]); 
%     axis(axis_range);
%     if ~isempty(fault_file)
%        for ii = 1:LS
%           slon = [lonf(ii) lonf(ii+LS)];
%           slat = [latf(ii) latf(ii+LS)];
% %           [xx,yy] = utm2ll(slon,slat,0,1);
%           [xx,yy] = ll2xy(slon,slat,ref_lon);
%           xs = (xx - xo) ./ 1000;
%           ys = (yy - yo) ./ 1000;
%           line(xs,ys,'color','black','linewidth',1.5);
%        end
%     end
%     
%     set(gcf,'PaperPositionMode','auto');
end