function plot_gps_model(gps_mat_file,gps_model,varargin)
   
   iter_step = 0;
   fault_file = [];
   lon_eq = -117.5; 
   lat_eq = 35.5;
   [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);   
   data_type = 'cont';   % default is cGPS
   site_file = [];
   model_type = 'okada';
   
   if ~isempty(varargin)
       for CC = 1:floor(length(varargin)/2)
           try
               switch lower(varargin{CC*2-1})
                   case 'iter_step'
                       iter_step = varargin{CC*2};
                   case 'fault'
                       fault_file = varargin{CC*2};
                   case 'data_type'
                       data_type = varargin{CC*2};
                   case 'site_name'
                       site_file = varargin{CC*2};
                   case 'model_type'
                       model_type = varargin{CC*2};
               end
           catch
               error('Unrecognized Keyword');
           end
       end
   end
   
   % read GPS data
   data = load(gps_mat_file);
   xin = data.data_gps(:,1) / 1000;
   yin = data.data_gps(:,2) / 1000;
   ux = data.data_gps(:,3);
   uy = data.data_gps(:,4);   
   
   if strcmp(data_type,'cont')
      sigx = data.data_gps(:,6);
      sigy = data.data_gps(:,7);
      Ngps = length(gps_model) / 3;
   elseif strcmp(data_type,'survey')
      sigx = data.data_gps(:,5);
      sigy = data.data_gps(:,6);
      Ngps = length(gps_model) / 2;
   else
      error('There is something wrong with the data type');
   end
   modelx = gps_model(1:Ngps);
   modely = gps_model(Ngps+1:2*Ngps);
   
   % read site name
   if ~isempty(site_file)
      fid = fopen(site_file);
      tmp = textscan(fid,'%s\n');
      fclose(fid);
      name_list = tmp{1}; 
      if length(name_list) ~= length(xin)
          error('Wrong number of GPS campaign sites');
      end
   end   
   
%    h0=figure('units','normalized','outerposition',[0 0 1 1]);  hold on
%    set(h0,'renderer','painters'); 
   h0 = figure; hold on
   quiver(xin,yin,modelx,modely,'Color','red','LineWidth',2,'AutoScale','off');
   quiver(xin,yin,ux,uy,'Color','blue','LineWidth',2,'AutoScale','off');
   % plot the error ellipse
   for kk = 1:length(xin)
       MU = double([xin(kk)+ux(kk),yin(kk)+uy(kk)]);
       plotEllipse(3*sigx(kk),3*sigy(kk),MU);
       if ~isempty(site_file)
          text(MU(1),MU(2),name_list{kk},'fontsize',15,'fontweight','bold','HorizontalAlignment','right');
       end
   end
   axis equal
   
   if strcmp(data_type,'cont')
      title('cGPS data fitting');
   elseif strcmp(data_type,'survey')
      title('Campaign GPS data fitting (Using campaign data)');
   else
      error('There is something wrong with the data type');
   end
   
   if ~isempty(fault_file)
       fault_trace = load(fault_file);
       lonf = [fault_trace(:,1);fault_trace(:,3)];  
       latf = [fault_trace(:,2);fault_trace(:,4)];
       LS = length(lonf) / 2;
       for ii = 1:LS
          slon = [lonf(ii) lonf(ii+LS)];
          slat = [latf(ii) latf(ii+LS)];
          [xx,yy] = utm2ll(slon,slat,0,1);
          xs = (xx - xo) ./ 1000;
          ys = (yy - yo) ./ 1000;
          line(xs,ys,'color','black','linewidth',1.5);
       end
   end   
   
   legend('Model','Data');
   set(gca,'Fontsize',20);   
   set(h0,'PaperPositionMode','auto');
   set(h0,'visible','off');
   saveas(h0,['GPS/',data_type,'_misfit_',num2str(iter_step)],'epsc');
   
   if strcmp(model_type,'okada')
       save_model_name = 'gps_model.mat';
   else
       save_model_name = 'gps_model_layered.mat';
   end   
   save(['GPS/',data_type,'_',save_model_name],'gps_model');
   
end