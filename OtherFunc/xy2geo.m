function [lat,lon]=xy2geo(lat0,lon0,xe,yn)
r=6371.100e3;
% r1d=(2.0*pi*r)/360.0;
% cov=pi/180.0;
% dlat=(yn/r1d)/cov;
% lat=lat0+dlat;
% c1d=0.5*2.0*pi*(r*cos(lat0*cov)+r*cos(lat*cov));
% dlon=(xe/c1d)/cov;
% lon=lon0+dlon;

cov=pi/180.0;
c=2*pi*r;
dlat=(yn/c)*360.0;
lat=lat0+dlat;
c1d=0.5*2.0*pi*(r*cos(lat0*cov)+r*cos(lat*cov));

dlon=(xe/c1d)*360.0;
lon=lon0+dlon;

