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

mask_near_field('dip_model/homogeneous/okada_best_model.mat','ASC64/branch_cut/three_subswath/los_samp3.mat','all_faults',150,'ASC64_los');
mask_near_field('dip_model/homogeneous/okada_best_model.mat','DES71/branch_cut/three_subswath/los_samp3.mat','all_faults',150,'DES71_los');
mask_near_field('dip_model/homogeneous/okada_best_model.mat','ASC64/offsets/los_samp3.mat','all_faults',150,'ASC64_rng');
mask_near_field('dip_model/homogeneous/okada_best_model.mat','DES71/offsets/los_samp3.mat','all_faults',150,'DES71_rng');
mask_near_field('dip_model/homogeneous/okada_best_model.mat','ALOS-2/T065/three_subswath/los_samp3.mat','all_faults',150,'T065_los');
mask_near_field('dip_model/homogeneous/okada_best_model.mat','ALOS-2/T066/three_subswath/los_samp3.mat','all_faults',150,'T066_los');
mask_near_field('dip_model/homogeneous/okada_best_model.mat','Cosmo_Skymed/ASC_offsets/los_samp3.mat','all_faults',150,'CSK_ASC');
mask_near_field('dip_model/homogeneous/okada_best_model.mat','Cosmo_Skymed/DES_offsets/los_samp3.mat','all_faults',150,'CSK_DES');
