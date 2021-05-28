function write_insar_sample(insar_data,file_out);
% print out the sampled insar data to the file
%
% Usage: write_insar_sample(insar_data,file_out);
% insar_data=[xout,yout,zout,ve,vn,vz,Npt,rms_out,xx1,xx2,yy1,yy2];
%
% by Kang Wang on 08/28/2015
%

format long

f_data_out=fopen(file_out,'w');
fprintf(f_data_out,'%15.4f %15.4f %13.4f %8.4f %8.4f %8.4f %8d %10.4f %12.4f %12.4f %12.4f %12.4f\n',insar_data');
fclose(f_data_out);

