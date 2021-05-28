clc
clear
format long

[x1,y1,z1]=grdread2('phase_patch_mask.grd');
[X,Y]=meshgrid(x1,y1);
X=double(X);
Y=double(Y);

indx_good=~isnan(z1);
xgood=double(X(indx_good));
ygood=double(Y(indx_good));
zgood=double(z1(indx_good));

zout=griddata(xgood,ygood,zgood,X,Y,'nearest');
zout(zout>pi)=pi;
zout(zout<-pi)=-pi;

grdwrite2(x1,y1,zout,'phase_patch_interp.grd');
quit;
