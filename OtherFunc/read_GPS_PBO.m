function [date_gps,n,e,u,sig_n,sig_e,sig_u]=read_GPS_PBO(fname);
% 
% read GPS time seris data in PBO format
%
%
% Usage:
% [date_gps,n,e,u,sig_n,sig_e,sig_u]=read_GPS_PBO(fname);
%
% by Kang Wang in June. 2016

fid=fopen(fname,'r');
head = 1; numhead = 0;
while head == 1
    line = fgetl(fid);
    if line(1) == ' '
        head = 0;
    else
        numhead = numhead + 1;
    end
end
fclose(fid);
[date_gps,n,e,u,sig_n,sig_e,sig_u] = textread(fname,'%s %*d %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %*f %*f %*f%*s', ...
    'headerlines',numhead);



