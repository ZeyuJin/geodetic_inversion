function plot_insar_model_resampled(sampled_data_file,los_model,varargin)
   data = load(sampled_data_file);
   losin = data.sampled_insar_data(:,3);
   los_res = losin - los_model;
   xin = data.sampled_insar_data(:,1) / 1000;
   yin = data.sampled_insar_data(:,2) / 1000;
   look_angle = data.sampled_insar_data(:,4:6);
   
   [filepath,~,~] = fileparts(sampled_data_file);
   indx = find(filepath == '/');
   label_name = filepath;
   save_name = filepath;
   if ~isempty(indx)
       if length(indx) == 2
          label_name = filepath(1:indx(2)-1); 
       else
          label_name = filepath(1:indx(1)-1); 
       end
       save_name = filepath(1:indx(1)-1); 
   end
   
   defo_max = 120; % cmax(losin);
   defo_min = -120; % min(losin);
   res_max = 20;
   iter_step = 0;
   fault_file = [];
   lon_eq = -117.5; 
   lat_eq = 35.5;
   ref_lon = lon_eq;
   axis_range = [50 150 -45 55];
   model_type = 'okada';
   
   if ~isempty(varargin)
       for CC = 1:floor(length(varargin)/2)
           try
               switch lower(varargin{CC*2-1})
                   case 'misfit_range'
                       res_max = varargin{CC*2};
                   case 'defo_max'
                       defo_max = varargin{CC*2};
                       defo_min = -defo_max;
                   case 'iter_step'
                       iter_step = varargin{CC*2};
                   case 'fault'
                       fault_file = varargin{CC*2};
                   case 'axis_range'
                       axis_range = varargin{CC*2};
                       if length(axis_range) ~= 4
                           error('Something wrong with axis range');
                       end
                   case 'model_type'
                       model_type = varargin{CC*2};  % homogenous or layered
                   case 'ref_lon'
                       ref_lon = varargin{CC*2};
                   case 'lonc'
                       lon_eq = varargin{CC*2};
                   case 'latc'
                       lat_eq = varargin{CC*2};
               end
           catch
               error('Unrecognized Keyword');
           end
       end
   end
   
   if ~isempty(fault_file)
       fault_trace = load(fault_file);
       lonf = [fault_trace(:,1);fault_trace(:,3)];  
       latf = [fault_trace(:,2);fault_trace(:,4)];
       LS = length(lonf) / 2;
   end
   
%    [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);
   [xo,yo] = ll2xy(lon_eq,lat_eq,ref_lon);
   
   sz = 30;
%    h0=figure('units','normalized','outerposition',[0 0 1 1]);
%    set(h0,'renderer','painters');
   h0 = figure;
   
   subplot('Position',[0.04 0.55 0.42 0.42]); hold on
   scatter(xin,yin,sz,losin,'filled');
   read_plot_fault_segment('SKFS_fault.txt','m');
   colormap jet
   colorbar
   title(['Sampled Data (',label_name,')']);
   set(gca,'Fontsize',20);
   caxis([defo_min defo_max]);
   axis(axis_range);
%    axis equal
   % plot the fault segments
   if ~isempty(fault_file)
       for ii = 1:LS
          slon = [lonf(ii) lonf(ii+LS)];
          slat = [latf(ii) latf(ii+LS)];
%           [xx,yy] = utm2ll(slon,slat,0,1);
          [xx,yy] = ll2xy(slon,slat,ref_lon);
          xs = (xx - xo) ./ 1000;
          ys = (yy - yo) ./ 1000;
          line(xs,ys,'color','black','linewidth',3);
       end
   end

   subplot('Position',[0.54 0.55 0.42 0.42]); hold on
   scatter(xin,yin,sz,los_model,'filled');
   read_plot_fault_segment('SKFS_fault.txt','m');
   colormap jet
   colorbar
   title('Model');
   set(gca,'Fontsize',20);
   caxis([defo_min defo_max]);
   axis(axis_range);
%    axis equal
   if ~isempty(fault_file)
       for ii = 1:LS
          slon = [lonf(ii) lonf(ii+LS)];
          slat = [latf(ii) latf(ii+LS)];
%           [xx,yy] = utm2ll(slon,slat,0,1);
          [xx,yy] = ll2xy(slon,slat,ref_lon);
          xs = (xx - xo) ./ 1000;
          ys = (yy - yo) ./ 1000;
          line(xs,ys,'color','black','linewidth',3);
       end
   end

   subplot('Position',[0.26 0.03 0.42 0.42]); hold on
   scatter(xin,yin,sz,los_res,'filled');
   read_plot_fault_segment('SKFS_fault.txt','m');
   colormap jet
   colorbar
   title('Residual');
   set(gca,'Fontsize',20);
   caxis([-res_max res_max]);       % center with zero
   axis(axis_range);
%    axis equal
   if ~isempty(fault_file)
       for ii = 1:LS
          slon = [lonf(ii) lonf(ii+LS)];
          slat = [latf(ii) latf(ii+LS)];
%           [xx,yy] = utm2ll(slon,slat,0,1);
          [xx,yy] = ll2xy(slon,slat,ref_lon);
          xs = (xx - xo) ./ 1000;
          ys = (yy - yo) ./ 1000;
          line(xs,ys,'color','black','linewidth',3);
       end
   end
   
   % save the sampled model for future use
   sampled_model = double([xin*1000,yin*1000,los_model,look_angle]);   % same format with sampled data
   if strcmp(model_type,'okada')
       save_model_name = 'los_model.mat';
   else
       save_model_name = 'los_model_layered.mat';
   end
   save([filepath,'/',save_model_name],'sampled_model');
   
   set(h0,'PaperPositionMode','auto');
%    set(h0,'visible','off');
%    saveas(h0,[filepath,'/',save_name,'_misfit_',num2str(iter_step)],'epsc');
   
%    % save the residual
%    if iter_step == 3
%        residual = [xin,yin,los_res];
%        save([filepath,'/',save_name,'_residual.mat'],'residual');
%    end
end