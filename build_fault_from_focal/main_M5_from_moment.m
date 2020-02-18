close all
clc
clear
addpath(genpath('/home/zeyu/software/codes_utilities'));

% focal mechanism information
% lon_eq = -117.575;   % origin of hypocenter
% lat_eq = 35.760;

lon_eq = -117.566916;  % center of aftershocks
lat_eq = 35.768305;
depth = 7e3;
strike = 218.68;
dip = 68;
Mw = 5.4;
sigma_w = 2e3;   % width of gaussian function
sigma_l = 2e3;
cmax = 25;       % colorbar range (25cm)

build_fault_M5(lon_eq,lat_eq,depth,strike,dip,Mw, ...
        'sig_len',sigma_l,'sig_width',sigma_w,'color_range',cmax);
         