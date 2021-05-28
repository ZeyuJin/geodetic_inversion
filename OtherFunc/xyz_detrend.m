function [zout,a,b,c]=xyz_detrend(xin,yin,zin);
% estimate and remove a trend from selected points
% Usage: [zout,a,b,c]=xyz_detrend(xpt,ypt,zpt);
%
% ztrend = a*xin+b*yin+c;
% zout  ------ zin - ztrend


indx_good=~isnan(zin);
x_good=xin(indx_good);
y_good=yin(indx_good);
z_good=zin(indx_good);

A=[x_good,y_good,ones(size(x_good))];
m=pinv(A)*z_good;
zfit=xin*m(1)+yin*m(2)+m(3);
zout=zin-zfit;

a=m(1);
b=m(2);
c=m(3);
