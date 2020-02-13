% close all
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

% %% inversion based on the best fault geometry
% dip_change_id = [1 2 3 4 5 7];
% H_offset = [-0.6 -0.2 0.5 0.5 0.6 -0.4] * 1e3;
% dip_angle = acosd(H_offset ./ 10e3);
% 
% slip_model_vs = load_fault_one_plane('fault_geometry/model_A_plane','dip_change_id',dip_change_id,'dip',dip_angle);
% %% re-add the top three layers of eastern branches
% % ID = 3;  N_layer_cut = 3;
% % [~,nL] = smoo1_each_plane(slip_model_vs);
% % patch_this_ID = nL(ID,:);
% % total_patch_before = sum(sum(nL(1:ID-1,:),2));
% % patch_3_layer_this_ID = 1:sum(patch_this_ID(1:N_layer_cut));
% % patch_disappear = total_patch_before + patch_3_layer_this_ID;
% % patch_rest_layer_this_ID = max(patch_disappear)+1:sum(sum(nL(1:ID,:),2));
% % % recompute the patch number and layer
% % slip_model_vs(patch_rest_layer_this_ID,2) = slip_model_vs(patch_rest_layer_this_ID,2) - length(patch_disappear);
% % slip_model_vs(patch_rest_layer_this_ID,3) = slip_model_vs(patch_rest_layer_this_ID,3) - N_layer_cut;
% % slip_model_vs(patch_disappear,:) = [];
% 
% curr_id = max(slip_model_vs(:,1));
% slip_model_ds = load_fault_dip_shallow('fault_geometry/model_A_splay','fault_id',curr_id,'dip_change_id',[1:7],'dip_depth',[3.2e3,3.2e3,2e3,4.2e3,4.2e3,4.2e3,3e3]);
% % slip_model = [slip_model_vs;slip_model_ds];
% % slip_model(:,2)=[1:size(slip_model,1)]';
% % show_slip_model(slip_model);
% % 
% % % nflt = max(slip_model(:,1));
% % % tSm = zeros(1,nflt+1);
% % % fault_id = slip_model(:,1);
% % % for i=1:nflt
% % %     tSm(i+1) = length(find(fault_id == i));
% % % end
% % 
% % % final inversion with resampled data
fault_file = 'fault_geometry/modelA_faults';
% segment_file = 'fault_geometry/model_A_seg_connect';     intersect_file = 'fault_geometry/model_A_dip_connect';
% iter_step = '3_detrend';
% model_type = 'okada';
% [slip_model,rms,model_roughness] = make_fault_from_insar(slip_model_vs,slip_model_ds,iter_step,'model_type',model_type, ...
%                   'fault',fault_file,'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'shallow_dip_id',[9:13]); 
% % save('dip_model/homogeneous/okada_best_model.mat','slip_model');
% % save('dip_model/layered/layered_edcmp.mat','slip_model');
              
%% test the data sensitivity
% dataset_sensitivity('ALOS2',slip_model_vs,slip_model_ds,iter_step, ...
%                     'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'shallow_dip_id',[9:13]);
% dataset_sensitivity('RNG',slip_model_vs,slip_model_ds,iter_step, ...
%                     'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'shallow_dip_id',[9:13]);
% dataset_sensitivity('AZO',slip_model_vs,slip_model_ds,iter_step, ...
%                     'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'shallow_dip_id',[9:13]);
% dataset_sensitivity('cgps',slip_model_vs,slip_model_ds,iter_step, ...
%                     'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'shallow_dip_id',[9:13]);
% dataset_sensitivity('camp_gps',slip_model_vs,slip_model_ds,iter_step, ...
%                     'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'shallow_dip_id',[9:13]);

% %% remove the ramp from the final model
% detrend_from_residual('ASC64/branch_cut/three_subswath','ramp_type','qu_ramp_7','dec',5,'fault',fault_file,'model_type','okada');
% detrend_from_residual('DES71/branch_cut/three_subswath','ramp_type','qu_ramp_7','dec',5,'fault',fault_file,'model_type','okada');
% detrend_from_residual('ALOS-2/T065/three_subswath','ramp_type','qu_ramp_5','dec',5,'fault',fault_file,'model_type','okada');
% detrend_from_residual('ALOS-2/T066/three_subswath','ramp_type','qu_ramp_7','dec',5,'fault',fault_file,'model_type','okada');
% detrend_from_residual('ASC64/offsets','ramp_type','noramp','dec',3,'fault',fault_file,'misfit_range',50,'model_type','okada');
% detrend_from_residual('DES71/offsets','ramp_type','bi_ramp','dec',3,'fault',fault_file,'misfit_range',50,'model_type','okada');
% detrend_from_residual('Cosmo_Skymed/ASC_offsets','ramp_type','bi_ramp','dec',2,'fault',fault_file,'misfit_range',50,'data_type','AZO','model_type','okada');
% detrend_from_residual('Cosmo_Skymed/DES_offsets','ramp_type','bi_ramp','dec',2,'fault',fault_file,'misfit_range',50,'data_type','AZO','model_type','okada');

% %% detrend the original and resampled dataset, and recompute the misfit
% detrend_final_dataset('ASC64/branch_cut/three_subswath','fault',fault_file,'axis_range',[-60,60,-40,80],'model_type','okada');
% detrend_final_dataset('DES71/branch_cut/three_subswath','fault',fault_file,'axis_range',[-60,60,-40,80],'model_type','okada');
% detrend_final_dataset('ALOS-2/T065/three_subswath','fault',fault_file,'axis_range',[-35,35,-5,60],'model_type','okada');
% detrend_final_dataset('ALOS-2/T066/three_subswath','fault',fault_file,'axis_range',[-45,35,-10,60],'model_type','okada');
% detrend_final_dataset('ASC64/offsets','fault',fault_file,'axis_range',[-30,20,0,50],'threshold',70,'misfit_range',50,'model_type','okada');
% detrend_final_dataset('DES71/offsets','fault',fault_file,'axis_range',[-30,20,5,50],'threshold',70,'misfit_range',50,'model_type','okada');
% detrend_final_dataset('Cosmo_Skymed/ASC_offsets','fault',fault_file,'axis_range',[-30,15,0,50],'threshold',70,'misfit_range',50,'model_type','okada');
% detrend_final_dataset('Cosmo_Skymed/DES_offsets','fault',fault_file,'axis_range',[-15,15,0,45],'threshold',70,'misfit_range',50,'model_type','okada');

plot_final_dataset('ASC64/branch_cut/three_subswath','axis_range',[-38 29 -5 55],'max_value',0.6001,'model_type','okada','title_name','Sentinel-1 (Asc 64) LOS','fault',fault_file);
plot_final_dataset('DES71/branch_cut/three_subswath','axis_range',[-38 29 -5 55],'max_value',0.6001,'model_type','okada','title_name','Sentinel-1 (Des 71) LOS','fault',fault_file,'label','def');
% plot_final_dataset('ALOS-2/T065/three_subswath','axis_range',[-38 29 -5 55],'max_value',0.6001,'model_type','okada','title_name','ALOS-2 (Asc 65) LOS','label','abc','fault',fault_file,'txt_pos',-1);
% plot_final_dataset('ALOS-2/T066/three_subswath','axis_range',[-38 29 -5 55],'max_value',0.6001,'model_type','okada','title_name','ALOS-2 (Asc 66) LOS','label','def','fault',fault_file,'txt_pos',-1);
% plot_final_dataset('ASC64/offsets','axis_range',[-38 29 -5 55],'threshold',70,'misfit_range',0.5,'max_value',1,'model_type','okada','title_name','Sentinel-1 (Asc 64) Range Offsets','label','abc','fault',fault_file,'txt_pos',-18);
% plot_final_dataset('DES71/offsets','axis_range',[-38 29 -5 55],'threshold',70,'misfit_range',0.5,'max_value',1,'model_type','okada','title_name','Sentinel-1 (Des 71) Range Offsets','label','def','fault',fault_file,'txt_pos',-18);
% plot_final_dataset('Cosmo_Skymed/ASC_offsets','axis_range',[-38 29 -5 55],'threshold',70,'misfit_range',0.5,'max_value',2,'model_type','okada','title_name','CSK (Asc 14) Azimuth Offsets','label','abc','fault',fault_file,'txt_pos',-13);
% plot_final_dataset('Cosmo_Skymed/DES_offsets','axis_range',[-38 29 -5 55],'threshold',70,'misfit_range',0.5,'max_value',2,'model_type','okada','title_name','CSK (Des 21) Azimuth Offsets','label','def','fault',fault_file,'txt_pos',-13);
