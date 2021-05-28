function GPS_ascii2mat(filepath,ascii_file,data_type)
% please provide the full path of ascii file of GPS data

   data = load([filepath,'/',ascii_file]);
   lon = data(:,1);
   lat = data(:,2);
   
   lon_eq = 72;
   lat_eq = 38.5;
   ref_lon = 71;
   [xo,yo] = ll2xy(lon_eq,lat_eq,ref_lon);
   [gps_x,gps_y] = ll2xy(lon,lat,ref_lon);
   
%    lon_eq = -117.5;
%    lat_eq = 35.5;
%    [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);   
%    [gps_x,gps_y] = utm2ll(lon,lat,0,1);
   gps_x = gps_x - xo;   gps_y = gps_y - yo;
   
   if strcmp(data_type,'cont')
       obs_three = data(:,3:5);
       sig_three = data(:,6:8);   
%        obs_three = obs_three .* 100;      % convert to cm
%        sig_three = sig_three .* 100;
       data_gps = double([gps_x,gps_y,obs_three,sig_three]);
       save([filepath,'/continuous_gps_3d.mat'],'data_gps');
       
   elseif strcmp(data_type,'survey')
       obs_two = data(:,3:4);
       sig_two = data(:,5:6);
%        obs_two = obs_two .* 100;
%        sig_two = sig_two .* 100;
       data_gps = double([gps_x,gps_y,obs_two,sig_two]);
       save([filepath,'/survey_gps_2d.mat'],'data_gps');
   else
       error('There is something wrong with the data type');
   end   
   
end