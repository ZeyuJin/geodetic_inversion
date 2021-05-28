function [X,Y]=ll2xy_grid(LON,LAT,lon0,lat0,lon_center);
% convert the matrix of longitude and latitude to cartesian coordinates
% 
% Usage: [X,Y]=ll2xy_grid(LAT,LAT,lon0,lat0,lon_center);
% 
% LON,LAT ----- matrix of longtitude and latitude
%
% X,Y    -----  matrix of cartesian coordinates
% 
% by Kang Wang in Nov. 2015
%

format long 
[xo,yo]=ll2xy(lon0,lat0,lon_center);

lon=LON(1,:);
lat=LAT(:,1);

x=lon;
y=lat;
[xutm1,yutm1]=ll2xy(x',y(1)*ones(size(x')),lon_center);
xe1=xutm1-xo;
yn1=yutm1-yo;


[xutm2,yutm2]=ll2xy(x(1)*ones(size(y)),y,lon_center);
xe2=xutm2-xo;
yn2=yutm2-yo;

[X,Y]=meshgrid(xe1,yn2);

