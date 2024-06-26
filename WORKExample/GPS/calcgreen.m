addpath(genpath('/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/Research/TurkeyEQ'))
load('GPS.mat', 'data_gps');
load('slip_model_optimized4.mat', 'slip_model');
G_raw = calc_green_gps_edcmp(slip_model,data_gps,'edgrnTK');
save('layered_green.mat', 'G_raw');