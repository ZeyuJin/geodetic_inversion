function [x,y]=geo2xy(lat,lon,lat0,lon0);

cov=pi/180.0;
r=6371000.0;

r1d=2*pi*r/360.0;

dlat=(lat-lat0);

dlon=(lon-lon0);
y=r1d*dlat;

lat_c=0.5*(lat0+lat)*cov;

r1d_lat=r1d*cos(lat_c);

x=r1d_lat*dlon;