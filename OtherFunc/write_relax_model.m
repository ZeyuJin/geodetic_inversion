function write_relax_model(relax_model,fname);
%
% Usage: write_relax_model(relax_model,fname);
%
% by Kang Wang on 10/27/2016

f_out=fopen(fname,'w');
fprintf(f_out,'%5d %12.5e %12.5f %12.5f %12.5f %12.5f %12.5f %8.2f %8.2f %8.2f\n',relax_model');
fclose(f_out);