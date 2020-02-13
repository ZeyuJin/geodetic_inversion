close all
clc
clear

data = load('/Users/zej011/Ridgecrest/data_resample/GPS/total_fixed.txt');
lon = data(:,1);
lat = data(:,2);

filename = 'cGPS_sites.kml';
kmlwritepoint(filename,lat,lon);