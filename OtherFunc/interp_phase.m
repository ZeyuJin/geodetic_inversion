clc
clear
format long

phase_in='phase_patch_mask.grd';
[x,y,z]=grdread2(phase_in);


indx_good=~isnan(z);
[X,Y]=meshgrid(x,y);
xgood=X(indx_good);
ygood=Y(indx_good);
zgood=z(indx_good);

zinterp=griddata(double(xgood),double(ygood),double(zgood),double(X),double(Y),'nearest');

%figure;
%h=imagesc(x,y,z);
%colorbar
%set(h,'alphadata',~isnan(z))

%figure;
%h=imagesc(x,y,zinterp);
%colorbar
%set(h,'alphadata',~isnan(zinterp))

file_out='phase_interp.grd';
grdwrite2(x,y,zinterp,file_out);
quit;
