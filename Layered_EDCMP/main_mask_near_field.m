close all
clc
clear

package_path = '/Users/zej011/work/zeyu/matlab/finite_fault_inversion';
addpath(genpath(package_path));
external_path = '/Users/zej011/work/zeyu/matlab/external_func';
addpath(genpath(external_path));
addpath('/Users/zej011/work/zeyu/matlab');
addpath('/Users/zej011/work/Kang_tutorial/candis');
addpath('/Users/zej011/work/Kang_tutorial/codes_utilities/matlab/igppsar');
setenv('PATH',[getenv('PATH'),':/Users/zej011/work/zeyu/cshell']);
setenv('PATH',[getenv('PATH'),':/opt/local/bin']);  % add the path of GMT

%% Ridgecrest earthquake
% mask_near_field('dip_model/homogeneous/okada_best_model_final.mat','ASC64/branch_cut/three_subswath/los_samp2_detrend.mat','all_faults',150,'ASC64_los');
% mask_near_field('dip_model/homogeneous/okada_best_model_final.mat','DES71/branch_cut/three_subswath/los_samp2_detrend.mat','all_faults',150,'DES71_los');
% mask_near_field('dip_model/homogeneous/okada_best_model_final.mat','ASC64/offsets/los_samp2_detrend.mat','all_faults',150,'ASC64_rng');
% mask_near_field('dip_model/homogeneous/okada_best_model_final.mat','DES71/offsets/los_samp2_detrend.mat','all_faults',150,'DES71_rng');
% mask_near_field('dip_model/homogeneous/okada_best_model_final.mat','ALOS-2/T065/three_subswath/los_samp2_detrend.mat','all_faults',150,'T065_los');
% mask_near_field('dip_model/homogeneous/okada_best_model_final.mat','ALOS-2/T066/three_subswath/los_samp2_detrend.mat','all_faults',150,'T066_los');
% mask_near_field('dip_model/homogeneous/okada_best_model_final.mat','Cosmo_Skymed/ASC_offsets/los_samp2_detrend.mat','all_faults',150,'CSK_ASC');
% mask_near_field('dip_model/homogeneous/okada_best_model_final.mat','Cosmo_Skymed/DES_offsets/los_samp2_detrend.mat','all_faults',150,'CSK_DES');
% mask_near_field('M1995_zeyu_homo.mat','ERS170/los_samp3.mat','plane_1995_fault',150,'ERS170');
% mask_near_field('M1995_zeyu_homo.mat','ERS442/los_samp3.mat','plane_1995_fault',150,'ERS442');
% mask_near_field('M1995_zeyu_homo.mat','ERS120/los_samp3.mat','plane_1995_fault',150,'ERS120');

%% Pamir earthquake
% mask_near_field('resample/improve_unwrap/homo_all_better_data.mat','ASC100/LOS2/los_samp1.mat','all_faults',200,'ASC100_los');
% mask_near_field('resample/improve_unwrap/homo_all_better_data.mat','DES5/LOS2/los_samp1.mat','all_faults',200,'DES5_los');
% mask_near_field('resample/improve_unwrap/homo_all_better_data.mat','ALOS2_SCAN2/los_samp1.mat','all_faults',200,'DES57_los');
% mask_near_field('resample/improve_unwrap/homo_all_better_data.mat','ALOS2_stripe/LOS/los_samp1.mat','all_faults',200,'ASC163_los');
mask_near_field('resample/improve_unwrap/homo_all_better_data.mat','ALOS2_stripe/MAI2/los_samp1.mat','all_faults',200,'ASC163_MAI');