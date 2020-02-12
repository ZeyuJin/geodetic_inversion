function plot_final_dataset(this_track,varargin)
% plot the resampled data, resampled model and full residual
    % default values
    lon_eq = -117.5;
    lat_eq = 35.5;
    res_max = 20/100;
    axis_range = [-100 100 -100 100];
    threshold = 50;
    max_value = 60/100;
    label = 'abc';
    model_type = [];
    title_name = [];
    fault_file = '';
    txt_left = -5;
    
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
                   case 'axis_range'
                       axis_range = varargin{CC*2};
                       if length(axis_range) ~= 4
                           error('Something wrong with axis range');
                       end
                   case 'threshold'
                       threshold = varargin{CC*2};
                   case 'max_value'
                       max_value = varargin{CC*2};
                   case 'label'
                       label = varargin{CC*2};
                       if length(label) ~= 3
                           error('Something wrong with labels');
                       end
                   case 'fault'
                       fault_file = varargin{CC*2};
                   case 'model_type'
                       model_type = varargin{CC*2};
                   case 'title_name'
                       title_name = varargin{CC*2};
                   case 'txt_pos'
                       txt_left = varargin{CC*2};
                       
               end
           catch
               error('Unrecognized Keyword');
           end
       end
    end 
   
   [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);
   % read fault data
   if ~isempty(fault_file)
      fault_trace = load(fault_file);
      lonf = [fault_trace(:,1);fault_trace(:,3)];  
      latf = [fault_trace(:,2);fault_trace(:,4)];
      LS = length(lonf) / 2;
   end
   
    % load resampled data and model
    if strcmp(model_type,'okada')
        save_data_name = 'los_samp3_detrend.mat';
        model_name = 'los_model.mat';
    else
        save_data_name = ['los_samp3_detrend_',model_type,'.mat'];
        model_name = 'los_model_layered.mat';
    end 
    data = load([this_track,'/',save_data_name]);
    los_samp = data.sampled_insar_data(:,3);
    XS = data.sampled_insar_data(:,1) / 1000;
    YS = data.sampled_insar_data(:,2) / 1000;   
    
    model = load([this_track,'/',model_name]);
    los_model = model.sampled_model(:,3);
    los_resid = los_samp - los_model;
    
    % delete some points which sampled on the boundary
    indx = find(los_resid > threshold);
    los_samp(indx) = NaN;
    los_model(indx) = NaN;
    
    [xfull,yfull,los_resid_full] = grdread2([this_track,'/los_residual_',model_type,'.grd']);
    [xmf,ymf] = meshgrid(xfull,yfull);
    [XU,YU] = utm2ll(xmf(:),ymf(:),0,1);
    XF = (XU - xo) ./ 1000;
    YF = (YU - yo) ./ 1000;
    XF = reshape(XF,size(xmf));
    YF = reshape(YF,size(ymf));
    los_resid_full(los_resid_full > threshold) = NaN;

    % plot the resampled data, model and full residual
    sz = 20;
    figure;
    subplot('Position',[0.027 0.5 0.323 0.4]); hold on
    h1 = scatter(XS,YS,sz,los_samp/100,'filled');
    shading flat
    colormap jet
%     colorbar('southoutside');
%     title('Sampled data (cm)');
%     text_x = axis_range(2)+12;
%     text_y = axis_range(4)-2;
    text(-37.5,60,['(',label(1),')'],'fontsize',20,'FontWeight','bold');
    set(gca,'Fontsize',20); 
    xlim(axis_range(1:2));
    ylim(axis_range(3:4));
    set(h1, 'YLimInclude', 'off');
    axis equal
    caxis([-max_value max_value]);
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
    set(gca,'box','on');
    
    % -18 for S1A-RNG
    % -5 for S1A-LOS
    % -1 for ALOS-2 LOS
    % -13 for CSK AZO
    text(txt_left,56,title_name,'fontsize',20,'FontWeight','bold');   
    
    subplot('Position',[0.35 0.5 0.323 0.4]); hold on
    h2 = scatter(XS,YS,sz,los_model/100,'filled');
    shading flat
    colormap jet
    cb1 = colorbar('southoutside');
    set(cb1,'position',[0.19,0.46,0.32,0.015]);
%     title(cb1,'cm');
%     lbpos = get(cb1,'title');
%     pos = get(lbpos,'position');
%     set(lbpos,'position',pos);
    set(gca,'Fontsize',20);
    xlim(axis_range(1:2));
    ylim(axis_range(3:4));
    set(h2, 'YLimInclude', 'off');
    text(-37.5,60,['(',label(2),')'],'fontsize',20,'FontWeight','bold');
    axis equal
    caxis([-max_value max_value]);
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
    set(gca,'YTick',[]);
    set(gca,'box','on');
    
    subplot('Position',[0.673 0.5 0.323 0.4]); hold on
    h3 = pcolor(XF,YF,los_resid_full/100);
    shading flat
    colormap jet
    cb2 = colorbar('southoutside');
    set(cb2,'position',[0.671,0.46,0.32,0.015]);
%     title('Full residual (cm)');
    set(gca,'Fontsize',20);
    text(-37.5,60,['(',label(3),')'],'fontsize',20,'FontWeight','bold');
    caxis([-res_max res_max]); 
    xlim(axis_range(1:2));
    ylim(axis_range(3:4));
    set(h3, 'YLimInclude', 'off');
    axis equal
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
    set(gca,'YTick',[]);
    set(gca,'box','on');
    read_plot_fault_segment('verified_rupture_trace.txt');
    read_plot_fault_segment('not_verified_rupture_trace.txt');
    
    set(gcf,'PaperPositionMode','auto');

end