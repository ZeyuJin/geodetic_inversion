function resamp_insar_data(slip_model_in, data_list,Nmin,Nmax,iter_step,varargin)
% resample the insar data with model predictions 
% method in Kang Wang and Yuri Fialko, GRL 2015
% written by Zeyu Jin on Sept. 2019
   set(0,'defaultAxesFontSize',15);
   
   % data_list defines each Interferogram track with such a format
   % the number means the number of data points to be resampled
   % ASC99    1200
   % DES106   1500
   
   % slip_model_in is the previous inverted model
   
   %% DEFINE THE DEFAULT VALUES
   iint=iter_step;
   lon_eq = -117.5;
   lat_eq = 35.5;
   ref_lon = lon_eq;
   fault_file = '';
   Nlook = 1;
   data_type = 'insar';
   
   %% read varargin values and assembly
   if ~isempty(varargin)
       for CC = 1:floor(length(varargin)/2)
           try
               switch lower(varargin{CC*2-1})
                   case 'lonc'
                       lon_eq = varargin{CC*2};
                   case 'latc'
                       lat_eq = varargin{CC*2};
                   case 'fault'
                       fault_file = varargin{CC*2};
                   case 'dec'
                       Nlook = varargin{CC*2};
                   case 'data_type'
                       data_type = varargin{CC*2};
                   case 'ref_lon'
                       ref_lon = varargin{CC*2};
               end
           catch
               error('Unrecognized Keyword\n');
           end
       end
   end
   
   %% to find how many tracks of data
   fid = fopen(data_list);
   tmp_txt = fgetl(fid);
   ntrack = 0;
   while tmp_txt ~= -1
       ntrack = ntrack + 1;
       tmp_txt = fgetl(fid);
   end
   fclose(fid);
   disp(['There are ',num2str(ntrack),' tracks of data using quadtree sampling strategy.']);
   
   %%  read txt file again to find those tracks and specify each sample regions
   track = cell(ntrack,1);   npt = zeros(ntrack,1);
   fid = fopen(data_list);
   tmp_txt = fgetl(fid);
   count = 0;
   while tmp_txt ~= -1
       count = count + 1;
       strs = strsplit(tmp_txt);
       track(count) = cellstr(strs{1});
       disp(count);
       disp(npt);
       disp(str2double(strs{2}));
       npt(count) = str2double(strs{2});  
       tmp_txt = fgetl(fid);
   end     
   fclose(fid);
   [xo,yo] = ll2xy(lon_eq,lat_eq,ref_lon);
   
   % % data covariance matrix
   % varigram=load('insar_varigram.mat');
   % sigma=varigram.sigma;
   % L=varigram.L;
   
   % iterative sample the data using model predictions
   for k=1:ntrack
       this_track=track{k};
       disp(['working on ',this_track]);
       this_npt=npt(k);
       
       % [x1,y1,demin]=grdread2([this_track,'/','dem_low.grd']);
       [x1,y1,losin]=grdread2([this_track,'/','los_ll_low','.grd']);   % in the unit of cm
       [x1,y1,ze]=grdread2([this_track,'/','look_e_low','.grd']);
       [x1,y1,zn]=grdread2([this_track,'/','look_n_low','.grd']);
       [x1,y1,zu]=grdread2([this_track,'/','look_u_low','.grd']);
       
       % multi-look to reduce the computation time
       if Nlook > 1
          % [lon1,lat1,deml] = multi_look(x1,y1,demin,Nlook,Nlook);
          [lon1,lat1,losl] = multi_look(x1,y1,losin,Nlook,Nlook);
          [lon1,lat1,zel] = multi_look(x1,y1,ze,Nlook,Nlook);
          [lon1,lat1,znl] = multi_look(x1,y1,zn,Nlook,Nlook);
          [lon1,lat1,zul] = multi_look(x1,y1,zu,Nlook,Nlook);
       else
           lon1 = x1;
           lat1 = y1;
           deml = demin;
           losl = losin;
           zel = ze;
           znl = zn;
           zul = zu;
       end
       
       [xm1,ym1] = meshgrid(lon1,lat1);
       [xutm,yutm] = ll2xy(xm1(:),ym1(:),ref_lon);
       xin = xutm - xo;
       yin = yutm - yo;
       xin = reshape(xin,size(xm1));
       yin = reshape(yin,size(ym1));
       
       % % put the slip model as an input parameter
%        slip_model_in = load('fault_M7.slip');
%        slip_model_in = load('fault_M6.slip');
%        tmp = load('resample/4_segments/homo_4data.mat');
%        tmp = load('resample/misfit_include_MAI/homo_better_data.mat');
%        slip_model_in = tmp.slip_model;
       
       if strcmp(data_type,'insar')
          los_model = slip2insar_okada(xin,yin,losl,zel,znl,zul,slip_model_in);   % fix the bug using multi-looked looking angles
       else
          los_model = slip2AZO_okada(xin,yin,losl,zel,znl,zul,slip_model_in);     % add module to compute AZO data
       end
       
       % resample the data using slip model
       % a bit similar to EM algorithm to make it diverge (E-step??  M-step?? Currently I don't know...)
       disp(size(lon1));
       disp(size(lat1));
       disp(size(los_model));
       [lon_model,lat_model,zout_model,Npt,rms_out,xx1,xx2,yy1,yy2]=make_insar_downsample(lon1,lat1,los_model,this_npt,Nmin,Nmax,'mean'); % same with downsample

       [lon_pt,lat_pt,zout]=make_look_downsample(lon1,lat1,losl,lon_model,lat_model,xx1,xx2,yy1,yy2);
       % [lon_pt,lat_pt,dem_out]=make_look_downsample(lon1,lat1,deml,lon_model,lat_model,xx1,xx2,yy1,yy2);
       [lon_pt,lat_pt,ve]=make_look_downsample(lon1,lat1,zel,lon_model,lat_model,xx1,xx2,yy1,yy2);
       [lon_pt,lat_pt,vn]=make_look_downsample(lon1,lat1,znl,lon_model,lat_model,xx1,xx2,yy1,yy2);
       [lon_pt,lat_pt,vz]=make_look_downsample(lon1,lat1,zul,lon_model,lat_model,xx1,xx2,yy1,yy2);
       
       % convert lon/lat to x/y
       [xutm,yutm]=ll2xy(lon_pt,lat_pt,ref_lon);
       xpt=xutm-xo;
       ypt=yutm-yo;
       
       indx_good=~isnan(zout);
       xpt=xpt(indx_good);
       ypt=ypt(indx_good);
       zout=zout(indx_good);
       % dem_out=dem_out(indx_good);
       ve=ve(indx_good);
       vn=vn(indx_good);
       vz=vz(indx_good);
       xx1=xx1(indx_good);
       yy1=yy1(indx_good);
       xx2=xx2(indx_good);
       yy2=yy2(indx_good);
           
       sampled_insar_data = double([xpt,ypt,zout,ve,vn,vz]);
   %     this_sig=sigma(k);
   %     this_L=L(k);
   %     covd = calc_insar_cov(xpt,ypt,this_sig,this_L); 
   %     save([this_track,'/','los_samp',num2str(iint),'.mat'],'insar_data','covd');
       save([this_track,'/los_samp',num2str(iint),'.mat'],'sampled_insar_data','rms_out');
       [hf,h1,h2]=plot_insar_sample_new(x1,y1,losin,zout,xx1,xx2,yy1,yy2,'fault',fault_file);
       set(hf,'PaperPositionMode','auto');
%        set(hf,'visible','off');
%       saveas(hf,[this_track,'/','los_samp',num2str(iint)],'epsc');
   end
   
end