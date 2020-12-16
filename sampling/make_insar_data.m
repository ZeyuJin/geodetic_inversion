function make_insar_data(data_list,Nmin,Nmax,varargin)
% downsample the insar data
% the data_list should be in the form: 
% path to data_file : number of points sampled : region size(optional)
% the script could apply uniform or quad-tree sampling based on curvature
% the varargin options:
% XY origin of UTM: lon_eq/lat_eq
% the option of sampling stretagy: uniform / quadtree?
set(0,'defaultAxesFontSize',15);

%% DEFINE THE DEFAULT VALUES
iint=num2str(0);          % number of step in iterative sampling
lon_eq = -117.5;
lat_eq = 35.5;
ref_lon = lon_eq;
sample_area = [-118.8 -116 34.3 36.9];       % for uniform sample in larger area
METHOD = 'quadtree';
fault_file = '';

%% read varargin values and assembly
if ~isempty(varargin)
    for CC = 1:floor(length(varargin)/2)
        try
            switch lower(varargin{CC*2-1})
                case 'method'
                    METHOD = varargin{CC*2};
                case 'lonc'
                    lon_eq = varargin{CC*2};
                case 'latc'
                    lat_eq = varargin{CC*2};
                case 'fault'
                    fault_file = varargin{CC*2};
                case 'area'
                    sample_area = varargin{CC*2};
                    if length(sample_area) ~= 4
                        error('There is something wrong with the input area!');
                    end
                case 'ref_lon'
                    ref_lon = varargin{CC*2};
            end
        catch
            error('Unrecognized Keyword');
        end
    end
end
if strcmp(METHOD,'quadtree')
    grd_file = 'los_clean_detrend.grd';
    file_suffix = 'low';
else
    grd_file = 'unwrap_clean_sample.grd';
    Nmax = Nmin;       % to force uniform sampling
    iint = '_uniform';
    file_suffix = 'uniform';
end

%% to find how many tracks of data
fid = fopen(data_list);
tmp_txt = fgetl(fid);
ntrack = 0;
while tmp_txt ~= -1
    ntrack = ntrack + 1;
    tmp_txt = fgetl(fid);
end
disp(['There are ',num2str(ntrack),' tracks of data using ',METHOD,' sampling strategy.']);
fclose(fid);

%%  read txt file again to find those tracks and specify each sample regions
track = cell(ntrack,1);   npt = zeros(ntrack,1);  region = zeros(ntrack,4); 
for ii = 1:ntrack
    region(ii,:) = sample_area;
end
fid = fopen(data_list);
tmp_txt = fgetl(fid);
count = 0;
while tmp_txt ~= -1
    count = count + 1;
    strs = strsplit(tmp_txt);
    track(count) = cellstr(strs{1});
    npt(count) = str2double(strs{2});
    num_of_strs = size(strs,2);
    % replaced by the region defined in the data_list
    if num_of_strs == 6, region(count,:) = str2num(char(strs{3:6})); end   
    tmp_txt = fgetl(fid);
end
fclose(fid);

% region(1,:) = [-118.8 -116 34.3 36.9];       % S1A/ASC64
% region(2,:) = [-118.8 -116 34.3 36.9];       % S1A/DES71
% region(3,:) = [-117.95 -116.60 35.2 36.3];   % ALOS-2/A065
% region(4,:) = [-118.3 -116.7 35.2 36.3];     % ALOS-2/A066
% region(5,:) = [-118.3 -116.7 35.2 36.3];     % ASC64/offsets
% region(6,:) = [-118.3 -116.7 35.2 36.3];     % DES71/offsets
% region(7,:) = [-118.3 -116.7 35.2 36.3];     % CSK/ASC/offsets
% region(8,:) = [-118.3 -116.7 35.2 36.3];     % CSK/DES/offsets

%% Do sampling
% sigma=zeros(ntrack,1);
% L=zeros(ntrack,1);
% [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);
[xo,yo] = ll2xy(lon_eq,lat_eq,ref_lon);
for k=1:ntrack
    this_track=track{k};
    [x1,y1,z1]=grdread2([this_track,'/',grd_file]);  % all data in CM unit
    los_this_track = z1;

%    los_noise=los_this_track_new(indx_out); %pixels used to estimate the covariance
%    [rr{k},vv{k},this_sig,this_L]=get_insar_varigram(xnoise,ynoise,double(los_noise),100,50e3,1000,2e6);
%    rnoise=rr{k};
%    vnoise=vv{k};
%    sigma(k)=this_sig;
%    L(k)=this_L;

%    xx=0:100:50e3;
%    yy=this_sig*exp(-xx/this_L);
%    hf=figure;
% %   set(hf,'visible','off')
%    plot(rnoise/1000,vnoise,'.','MarkerSize',10)
%    hold on
%    plot(xx/1000,yy,'r-','LineWidth',2);
%    xlabel('Distance (km)');
%    ylabel('Variance (cm^2)');
%    grid on
%    saveas(hf,[this_track,'/','los_variance'],'eps2c');
    
    this_npt=npt(k);
    
    xmin=region(k,1);
    xmax=region(k,2);
    ymin=region(k,3);
    ymax=region(k,4);
    
    indx_x=find(x1>=xmin & x1<=xmax);
    indx_y=find(y1>=ymin & y1<=ymax);
    xin=x1(indx_x);
    yin=y1(indx_y);
    losin=los_this_track(indx_y,indx_x);
    
    [xdem,ydem,zdem] = grdread2([this_track,'/dem_samp.grd']);
    demin = zdem(indx_y,indx_x);
    clear xdem ydem
    
    % make sure the grid size of looking angle same with phase (without multi-looking)
    [x2,y2,ze]=grdread2([this_track,'/look_e.grd']);
    [x3,y3,zn]=grdread2([this_track,'/look_n.grd']);
    [x4,y4,zu]=grdread2([this_track,'/look_u.grd']);
    ein=ze(indx_y,indx_x);
    nin=zn(indx_y,indx_x);
    uin=zu(indx_y,indx_x);
    clear x1 y1 ze zn zu
%
%     [lon_out,lat_out,los_out]=multi_look(x1,y1,los_this_track_new,5,5);
%     [lon_out,lat_out,ze_out]=multi_look(x2,y2,ze,5,5);    
%     [lon_out,lat_out,zn_out]=multi_look(x3,y3,zn,5,5);    
%     [lon_out,lat_out,zu_out]=multi_look(x4,y4,zu,5,5);    
%
%     grdwrite2(lon_out,lat_out,los_out,[this_track,'/los_ll_low.grd']);
%     grdwrite2(lon_out,lat_out,ze_out,[this_track,'/look_e_low.grd']);
%     grdwrite2(lon_out,lat_out,zn_out,[this_track,'/look_n_low.grd']);
%     grdwrite2(lon_out,lat_out,zu_out,[this_track,'/look_u_low.grd']);  
%  
   grdwrite2(xin,yin,losin,[this_track,'/los_ll_',file_suffix,'.grd']);
   grdwrite2(xin,yin,ein,[this_track,'/look_e_',file_suffix,'.grd']); 
   grdwrite2(xin,yin,nin,[this_track,'/look_n_',file_suffix,'.grd']);
   grdwrite2(xin,yin,uin,[this_track,'/look_u_',file_suffix,'.grd']);
   grdwrite2(xin,yin,demin,[this_track,'/dem_',file_suffix,'.grd']);

   % the resolution of geocoded insar data is about 100 meters
%    Nmin = 3;   Nmax = 150;  for quad-tree sampling
%    Nmin = Nmax = 15;        for uniform sampling
   [xout,yout,zout,Npt,rms_out,xx1,xx2,yy1,yy2]=make_insar_downsample(xin,yin,losin,this_npt,Nmin,Nmax,'mean');
%    [xutm_sar,yutm_sar] = utm2ll(xout,yout,0,1);
   [xutm_sar,yutm_sar] = ll2xy(xout,yout,ref_lon);
   xsar=xutm_sar-xo;
   ysar=yutm_sar-yo; 
%    covd = calc_insar_cov(xsar,ysar,this_sig,this_L);
%   
   [xe_out,ye_out,ve]=make_look_downsample(xin,yin,ein,xout,yout,xx1,xx2,yy1,yy2);
   [xn_out,yn_out,vn]=make_look_downsample(xin,yin,nin,xout,yout,xx1,xx2,yy1,yy2);
   [xu_out,yu_out,vz]=make_look_downsample(xin,yin,uin,xout,yout,xx1,xx2,yy1,yy2);
   [xdem_out,ydem_out,dem_out]=make_look_downsample(xin,yin,demin,xout,yout,xx1,xx2,yy1,yy2);
   
   sampled_insar_data=double([xsar,ysar,zout,ve,vn,vz]);     % save the downsampled insar data
   save([this_track,'/los_samp',iint,'.mat'],'sampled_insar_data','rms_out','dem_out');
   
   [h0,h1,h2]=plot_insar_sample_new(xin,yin,losin,zout,xx1,xx2,yy1,yy2,'fault',fault_file); 
   set(h0,'PaperPositionMode','auto');
%    set(h0,'visible','off');
   saveas(h0,[this_track,'/','los_samp',num2str(iint)],'epsc');
end