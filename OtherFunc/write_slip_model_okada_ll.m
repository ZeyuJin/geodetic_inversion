function write_slip_model_okada_ll(slip_model,slip_name);
% print out the slip model 
% 
% Usage: write_slip_model_okada_ll(slip_model,slip_name);
% 
% slip_model:  matrix containing the information of the slip model
% slip_name:   string of the file name 
% 

fid=fopen(slip_name,'w');
fprintf(fid,'%5d %5d %5d %15.5f %15.5f %10.2f %10.2f %10.2f %6.1f %6.1f %6.3f %8.3f %8.3f\n ',slip_model');
fclose(fid);
