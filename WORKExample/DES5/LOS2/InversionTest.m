%% Geodetic Inversion Simplified Codes
%11/2/2022, Xiaoyu Zou, x3zou@ucsd.edu
clear
clc
addpath('/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/geodetic_inversion-master/sign_mask')
addpath('/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/geodetic_inversion-master/WORK/ASC100/LOS2')
addpath('/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/geodetic_inversion-master/WORK/DES5/LOS2')
setenv('PATH',[getenv('PATH'),':/opt/homebrew/bin']);  % add the path of GMT


%% Step 1: data cleaning using clean_insar-data. Remove some near-field
%% unwrapping erros manually first.
clean_insar_data

%% Step 2: detrend the phase and remove the phase ambiguity
grdin = '/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/geodetic_inversion-master/WORK/DES5/LOS2/unwrap_clean_sample.grd' ;
grdout = '/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/geodetic_inversion-master/WORK/DES5/LOS2/los_clean_detrend.grd' ;

lonf = 73.215150;   latf = 37.754558;  %coordinate of pixel point
ref_lon = 71;   threshold = 0.5;
lonc = 72;  latc = 38.5; %coordinate of reference point (0,0)
remove_ref_from_grid(grdin,grdout,lonf,latf,ref_lon,threshold);

%% Step 3: apply quad-tree sampling to all detrended data 
data_list = 'detrend.list';     % LOS data
fault_file = 'all_faults';
Nmin = 2;    Nmax = 400;
make_insar_data(data_list,Nmin,Nmax,'method','quadtree','fault',fault_file,'area',[71.8 73.9 37.7 39.1],'ref_lon',71);


%% Step 4: Build the fault geometry
dip_change_id=1:5;%the array of fault ids that have dip angles not equal to 90 degrees
dip_angle=[87.7, 81.8, 85, 89.3, 89.3];%the array of fault ids that have dip angles not equal to 90 degrees
slip_model_vs = load_fault_one_plane(fault_file,'dip_change_id',dip_change_id,'dip',dip_angle,'lonc',72,'latc',38.5,'ref_lon',71,'len_top',1.2e3);
slip_model_ds=[];%set empty for most cases, because it's used to construct the geometry of "Y shape" or "flower structure" formed by shallow splay faults.

%% Step 5: inversion using first downsampled data
iter_step = 0;
segment_smooth_file = 'seg_connect';
intersect_smooth_file = [];
intersect_file=[];
segment_file = 'seg_connect';
[slip_model,~,~,~,~] = make_fault_from_insar(slip_model_vs,slip_model_ds,iter_step,'shallow_dip_id',[], ...
                  'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'fault',fault_file, ...
                  'lonc',72,'latc',38.5,'ref_lon',71,'model_type','okada');


%% Step 6: iterative sampling data using the model predictions (Wang and
%% Fialko, GRL 2015) 
iter_step = 1;  % usually just one iteration is enough to rule out samples on noisy pixels
los_list=data_list;
resamp_insar_data(slip_model,data_list,Nmin, Nmax, iter_step, 'fault', fault_file, 'dec',2, 'lonc',72, 'latc',38.5, 'ref_lon',71);

%% Step 7: inversion using resampled data
[slip_model,rms,model_roughness,insar_model1,insar_model2] = make_fault_from_insar(slip_model_vs,slip_model_ds,iter_step, ...
                     'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'fault',fault_file, ...
                     'lonc',72,'latc',38.5,'ref_lon',71);
compute_moment(slip_model,model_type);

%% Step 7.5: plotting the finite fault inversion and resampled data fitting
model_type = 'okada';
lon_eq = -117.5;
lat_eq = 35.5;
show_slip_model(slip_model,'misfit_range',400,'ref_lon',ref_lon,'lonc',lon_eq,'latc',lat_eq,'axis_range',[60 110 -45 25 -25 0]);
iint=iter_step;
plot_insar_model_resampled(['/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/geodetic_inversion-master/WORK/ASC100/LOS2/los_samp',num2str(iint),'.mat'],insar_model1,'iter_step',iint,'fault',fault_file,'model_type',model_type,'misfit_range',30,'defo_max',120,'ref_lon',ref_lon,'lonc',lon_eq,'latc',lat_eq);
plot_insar_model_resampled(['/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/geodetic_inversion-master/WORK/DES5/LOS2/los_samp',num2str(iint),'.mat'],insar_model2,'iter_step',iint,'fault',fault_file,'model_type',model_type,'misfit_range',30,'defo_max',120,'ref_lon',ref_lon,'lonc',lon_eq,'latc',lat_eq);
