function write_slip_model_fem(slip_model,slip_name);
% print the slip model (Node-based) to a file
%
%  Usage: write_slip_model_fem(slip_model,slip_name);
%  
%  slip_model: matrix containing the information of the slip model
%  slip_name: string of the file name
%
% by Kang Wang on 08/24/2015

fid=fopen(slip_name,'w');
fprintf(fid,'%5d %10.2f %10.2f %10.2f %8.3f %8.3f %8.3f %10.3f %10.3f\n',slip_model');
fclose(fid);


