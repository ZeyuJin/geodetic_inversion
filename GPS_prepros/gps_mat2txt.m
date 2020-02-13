function gps_mat2txt(GPS_data_file,GPS_model_file,dim,write_txt)
   d1 = load(GPS_data_file);
   lon = d1(:,1);  lat = d1(:,2);
   
   d2 = load(GPS_model_file);
   gps_model = d2.gps_model;
   
   Ngps = length(gps_model) / dim;
   Ue = gps_model(1:Ngps);
   Un = gps_model(Ngps+1:2*Ngps);
   error = zeros(Ngps,3);
   if dim == 3
      data = [lon,lat,Ue,Un,error]';
   else
      error(:,3) = 0.001;
      data = [lon,lat,Ue,Un,error]';
   end
   fid = fopen(write_txt,'w');
   fprintf(fid,'%f %f %f %f %f %f %f\n',data);
   fclose(fid);
end