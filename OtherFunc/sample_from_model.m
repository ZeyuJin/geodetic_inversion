%% Sample the data from model
%% 1. Save InSAR model 
%% 2. Sample the original data based on InSAR model
%% 3. Now the outputs are 'los_samp1'
% Xiaoyu Zou, 03/15/2025

clear
addpath(genpath('/Volumes/T7/Research/Tingri_project/geodetic_inversion-master'));
fault_file={'./Google_Earth_Data/fault_main.txt','./Google_Earth_Data/fault_sub.txt'};
lonc=87.5;
latc=28.7;
ref_lon=87.5;
iint=1;
data_list='sample_list_model2.txt'; %Column: data directory 1; data file 2; data type 3; npt (not used now) 4; bounds 5-8; Nmin, Nmax (size in km): 9-10; variance threshold: 11
nan_frac_max=1; % default: 0.7
% make_insar_data_from_model(data_list,'method','quadtree','lonc',lonc,'latc',latc,'ref_lon',ref_lon,'fault',fault_file,'iint',1,'nan_frac_max',nan_frac_max)
%mesh_trace={'mesh_trace_main.txt','mesh_trace_sub.txt'};
trace1=load('mesh_trace_main.txt');
trace2=load('mesh_trace_sub.txt');
faults={trace1,trace2};

make_insar_data_from_model(data_list,'method','quadtree','lonc',lonc,...
    'latc',latc,'ref_lon',ref_lon,'fault',fault_file,...
    'iint',iint,'nan_frac_max',nan_frac_max,'low',1,...
    'mesh_trace',trace2,'sign_clean',false)
