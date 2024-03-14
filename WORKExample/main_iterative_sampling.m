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
setenv('PATH',[getenv('PATH'),':/usr/local/bin']);  % add the path of GMT

ref_lon = 71;

% % % apply quad-tree sampling to all detrended data (LOS + RNG)
% fault_file = 'all_faults';
% area = [71.8 73.9 37.7 39.1];
% los_list = 'detrend.list';
% N1min = 4;   N1max = 300;
% make_insar_data(los_list,N1min,N1max,'method','quadtree','fault',fault_file,'ref_lon',71,'area',area,'lonc',72,'latc',38.5);
% 
% rng_list = 'ASC100_rng';
% N2min = 4;   N2max = 400;
% make_insar_data(rng_list,N2min,N2max,'method','quadtree','fault',fault_file,'ref_lon',71,'area',area,'lonc',72,'latc',38.5);

% ALOS2_offsets = 'ALOS2_los';
% N3min = 3;   N3max = 300;
% make_insar_data(ALOS2_offsets,N3min,N3max,'method','quadtree','fault',fault_file,'ref_lon',71,'area',area,'lonc',72,'latc',38.5);
% 
% ALOS2_offsets = 'ALOS2_rng';
% N3min = 3;   N3max = 300;
% make_insar_data(ALOS2_offsets,N3min,N3max,'method','quadtree','fault',fault_file,'ref_lon',71,'area',area,'lonc',72,'latc',38.5);

% ALOS2_azi = 'ALOS2_azo';
% N3min = 6;   N3max = 500;
% make_insar_data(ALOS2_azi,N3min,N3max,'method','quadtree','fault',fault_file,'ref_lon',71,'area',area,'lonc',72,'latc',38.5);

% ALOS2_SCAN_test = 'ALOS2_SCAN_test';
% N3min = 3;   N3max = 300;
% make_insar_data(ALOS2_SCAN_test,N3min,N3max,'method','quadtree','fault',fault_file,'ref_lon',71,'area',area,'lonc',72,'latc',38.5);
% 
% ASC100_test = 'ASC100_test';
% N3min = 3;  N3max = 300;
% make_insar_data(ASC100_test,N3min,N3max,'method','quadtree','fault',fault_file,'ref_lon',71,'area',area,'lonc',72,'latc',38.5);

% DES5_test = 'DES5_test';
% N3min = 3;  N3max = 300;
% make_insar_data(ASC100_test,N3min,N3max,'method','quadtree','fault',fault_file,'ref_lon',71,'area',area,'lonc',72,'latc',38.5);

%% inversion based on the best fault geometry
dip_change_id = 1:5;

% first inversion based on the downsample data
iter_step = 1; 
segment_file = 'seg_connect';  
intersect_file = [];
slip_model_ds = [];

% for i = 1:length(dip_NE)
    
%     tmp = dip_NE(i);
%     tmp = 89.3;
%     dip_angle = [87.7, 81.8, 81.8, tmp, 89.3];

dip_angle = [87.7, 81.8, 85, 89.3, 89.3];
% dip_change_id = 1:3;
% dip_angle = [87.7, 81.8, 89.3];

    fault_file = 'all_faults';
% fault_file = 'fault_3_seg_east';
    slip_model_vs = load_fault_one_plane(fault_file,'dip_change_id',dip_change_id,'dip',dip_angle, ...
                                     'lonc',72,'latc',38.5,'ref_lon',71,'len_top',1.2e3);
                                 
% slip_model = [slip_model_vs;slip_model_ds];
% show_slip_model(slip_model,'ref_lon',ref_lon);
    [slip_model,~,~] = make_fault_from_insar(slip_model_vs,slip_model_ds,iter_step,'shallow_dip_id',[], ...
                  'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'fault',fault_file, ...
                  'lonc',72,'latc',38.5,'ref_lon',71,'model_type','okada');
% end
               
% % % iterative resampling and inversion
% fault_file = 'all_faults';
% % % N1min = 2;   N1max = 300;
% N2min = 3;   N2max = 300;
% % % los_list = 'detrend.list';   % rng_list = 'rng.list'; 
% % 
% % % ASC100_rng = 'ASC100_rng';
% % % DES5_rng = 'DES5_rng';
% % % ALOS2_rng = 'ALOS2_rng';      
% % % ALOS2_los = 'ALOS2_los';
% % % ALOS2_azo = 'ALOS2_azo';
% % ALOS2_SCAN_test = 'ALOS2_SCAN_test';
% % ASC100_test = 'ASC100_test';
% DES5_test = 'DES5_test';
% % segment_file = 'seg_connect';     intersect_file = [];
% %  
% for iter_step = 1:1     % first two resampling steps
% % SF = [1e-3 5e-3 1e-2 2e-2 3e-2 4e-2 5e-2 6e-2 7e-2 8e-2 9e-2 1e-1 2e-1 3e-1 4e-1 5e-1 1];   % smoothness test
% % RMS = zeros(size(SF));
% % Rough = zeros(size(SF));
% % for jj = 1:length(SF)
%     disp(['Iterative sampling at No.',num2str(iter_step),' step']);
%     resamp_insar_data(DES5_test,N2min,N2max,iter_step,'fault',fault_file,'dec',1,'lonc',72,'latc',38.5,'ref_lon',71);
%     
% %     resamp_insar_data(ASC100_test,N2min,N2max,iter_step,'fault',fault_file,'dec',1,'lonc',72,'latc',38.5,'ref_lon',71);
% %     resamp_insar_data(los_list,N1min,N1max,iter_step,'fault',fault_file,'dec',2,'lonc',72,'latc',38.5,'ref_lon',71);
% %     resamp_insar_data(rng_list,N2min,N2max,iter_step,'fault',fault_file,'dec',2,'lonc',72,'latc',38.5,'ref_lon',71);
% %     resamp_insar_data(ALOS2_los,N2min,N2max,iter_step,'fault',fault_file,'dec',1,'lonc',72,'latc',38.5,'ref_lon',71);
% %     resamp_insar_data(ALOS2_SCAN_test,N2min,N2max,iter_step,'fault',fault_file,'dec',1,'lonc',72,'latc',38.5,'ref_lon',71);
% %     resamp_insar_data(DES5_rng,N2min,N2max,iter_step,'fault',fault_file,'dec',2,'lonc',72,'latc',38.5,'ref_lon',71);
% %     resamp_insar_data(ALOS2_azo,N2min,N2max,iter_step,'fault',fault_file,'dec',1,'lonc',72,'latc',38.5,'ref_lon',71,'data_type','AZO');
% %     lambda = SF(jj);
% %     [slip_model,rms,model_roughness] = make_fault_from_insar(slip_model_vs,slip_model_ds,iter_step, ...
% %                     'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'fault',fault_file, ...
% %                     'lonc',72,'latc',38.5,'ref_lon',71);
% %     RMS(jj) = rms;
% %     Rough(jj) = model_roughness;
% % end
% % L_data = [RMS',Rough'];
% % save('L_curve.mat','L_data');
% end

% % final inversion with resampled data
% % fault_file = 'all_faults';
% % N1min = 3;   N1max = 300;
% % N2min = 3;   N2max = 400;
% % N3min = 3;   N3max = 400;
% % los_list = 'insar.list';   rng_list = 'RNG.list';   AZO_list = 'CSK.list';
% fault_file = 'fault_geometry/modelA_faults';
% segment_file = 'fault_geometry/model_A_seg_connect';     intersect_file = 'fault_geometry/model_A_dip_connect';
% 
% iter_step = '3_mask';
% % disp(['Iterative sampling at No.',num2str(iter_step),' step']);
% % resamp_insar_data(los_list,N1min,N1max,iter_step,'fault',fault_file,'dec',1);
% % resamp_insar_data(rng_list,N2min,N2max,iter_step,'fault',fault_file,'dec',1);
% % resamp_insar_data(AZO_list,N3min,N3max,iter_step,'fault',fault_file,'dec',1,'data_type','AZO');
% [slip_model,rms,model_roughness] = make_fault_from_insar(slip_model_vs,slip_model_ds,iter_step,'model_type','edcmp', ...
%                   'fault',fault_file,'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'shallow_dip_id',[9:13]); 
% % save('dip_model/layered/tmp.mat','slip_model');


%% remove the ramp from the final model and plot final results
% detrend_from_residual('ASC64/branch_cut/three_subswath','ramp_type','qu_ramp_7','dec',5,'fault',fault_file,'model_type','layered');
% detrend_from_residual('DES71/branch_cut/three_subswath','ramp_type','qu_ramp_7','dec',5,'fault',fault_file,'model_type','layered');
% detrend_from_residual('ALOS-2/T065/three_subswath','ramp_type','qu_ramp_5','dec',5,'fault',fault_file,'model_type','layered');
% detrend_from_residual('ALOS-2/T066/three_subswath','ramp_type','qu_ramp_7','dec',5,'fault',fault_file,'model_type','layered');
% detrend_from_residual('ASC64/offsets','ramp_type','noramp','dec',3,'fault',fault_file,'misfit_range',50,'model_type','layered');
% detrend_from_residual('DES71/offsets','ramp_type','bi_ramp','dec',3,'fault',fault_file,'misfit_range',50,'model_type','layered');
% detrend_from_residual('Cosmo_Skymed/ASC_offsets','ramp_type','bi_ramp','dec',2,'fault',fault_file,'misfit_range',50,'data_type','AZO','model_type','layered');
% detrend_from_residual('Cosmo_Skymed/DES_offsets','ramp_type','bi_ramp','dec',2,'fault',fault_file,'misfit_range',50,'data_type','AZO','model_type','layered');
% detrend_from_residual('ALOS2_stripe/MAI','ramp_type','bi_ramp','dec',1,'lonc',72,'latc',38.5,'ref_lon',71,'misfit_range',200,'fault',fault_file,'data_type','AZO');
% detrend_from_residual('ALOS2_SCAN','ramp_type','qu_ramp_5','dec',5,'fault',fault_file,'lonc',72,'latc',38.5,'ref_lon',71,'misfit_range',60);
% detrend_from_residual('ALOS2_SCAN2','ramp_type','qu_ramp_5','dec',5,'fault',fault_file,'lonc',72,'latc',38.5,'ref_lon',71,'misfit_range',60);
% detrend_from_residual('ASC100/LOS2','ramp_type','noramp','dec',5,'fault',fault_file,'lonc',72,'latc',38.5,'ref_lon',71,'misfit_range',100);
% detrend_from_residual('DES5/LOS2','ramp_type','noramp','dec',5,'fault',fault_file,'lonc',72,'latc',38.5,'ref_lon',71,'misfit_range',100);
% detrend_from_residual('ALOS2_stripe/MAI2','ramp_type','noramp','dec',2,'fault',fault_file,'lonc',72,'latc',38.5,'ref_lon',71,'misfit_range',100,'data_type','AZO');
% detrend_from_residual('ALOS2_stripe/LOS','ramp_type','noramp','dec',2,'fault',fault_file,'lonc',72,'latc',38.5,'ref_lon',71,'misfit_range',100);

% %% detrend the original and resampled dataset, and recompute the misfit
% detrend_final_dataset('ASC64/branch_cut/three_subswath','fault',fault_file,'axis_range',[-60,60,-40,80],'model_type','layered');
% detrend_final_dataset('DES71/branch_cut/three_subswath','fault',fault_file,'axis_range',[-60,60,-40,80],'model_type','layered');
% detrend_final_dataset('ALOS-2/T065/three_subswath','fault',fault_file,'axis_range',[-35,35,-5,60],'model_type','layered');
% detrend_final_dataset('ALOS-2/T066/three_subswath','fault',fault_file,'axis_range',[-45,35,-10,60],'model_type','layered');
% detrend_final_dataset('ASC64/offsets','fault',fault_file,'axis_range',[-30,20,0,50],'threshold',70,'misfit_range',50,'model_type','layered');
% detrend_final_dataset('DES71/offsets','fault',fault_file,'axis_range',[-30,20,5,50],'threshold',70,'misfit_range',50,'model_type','layered');
% detrend_final_dataset('Cosmo_Skymed/ASC_offsets','fault',fault_file,'axis_range',[-30,15,0,50],'threshold',70,'misfit_range',50,'model_type','layered');
% detrend_final_dataset('Cosmo_Skymed/DES_offsets','fault',fault_file,'axis_range',[-15,15,0,45],'threshold',70,'misfit_range',50,'model_type','layered');
% detrend_final_dataset('ALOS2_stripe/MAI','lonc',72,'latc',38.5,'ref_lon',71);
% detrend_final_dataset('ALOS2_SCAN','lonc',72,'latc',38.5,'ref_lon',71);

% plot_final_dataset('ASC64/branch_cut/three_subswath','axis_range',[-38 29 -5 55],'model_type','layered');
% plot_final_dataset('DES71/branch_cut/three_subswath','axis_range',[-38 29 -5 55],'model_type','layered');
% plot_final_dataset('ALOS-2/T065/three_subswath','axis_range',[-38 29 -5 55],'model_type','layered');
% plot_final_dataset('ALOS-2/T066/three_subswath','axis_range',[-38 29 -5 55],'model_type','layered');
% plot_final_dataset('ASC64/offsets','axis_range',[-38 29 -5 55],'threshold',70,'misfit_range',50,'max_value',100,'model_type','layered');
% plot_final_dataset('DES71/offsets','axis_range',[-38 29 -5 55],'threshold',70,'misfit_range',50,'max_value',100,'model_type','layered');
% plot_final_dataset('Cosmo_Skymed/ASC_offsets','axis_range',[-38 29 -5 55],'threshold',70,'misfit_range',50,'max_value',200,'model_type','layered');
% plot_final_dataset('Cosmo_Skymed/DES_offsets','axis_range',[-38 29 -5 55],'threshold',70,'misfit_range',50,'max_value',200,'model_type','layered');
% plot_final_dataset('ALOS2_SCAN','fault',fault_file,'lonc',72,'latc',38.5,'ref_lon',71,'misfit_range',30,'max_value',120,'threshold',60,'axis_range',[45 115 -45 25]);
% plot_final_dataset('ALOS2_SCAN2','fault',fault_file,'lonc',72,'latc',38.5,'ref_lon',71,'misfit_range',30,'max_value',120,'threshold',60,'axis_range',[45 115 -45 25]);
% plot_final_dataset('ASC100/LOS2','fault',fault_file,'lonc',72,'latc',38.5,'ref_lon',71,'misfit_range',30,'max_value',120,'threshold',60,'axis_range',[45 115 -45 25]);
% plot_final_dataset('DES5/LOS2','fault',fault_file,'lonc',72,'latc',38.5,'ref_lon',71,'misfit_range',30,'max_value',120,'threshold',60,'axis_range',[45 115 -45 25],'label','def');
% plot_final_dataset('ALOS2_stripe/MAI2','fault',fault_file,'lonc',72,'latc',38.5,'ref_lon',71,'misfit_range',60,'max_value',120,'threshold',60,'axis_range',[45 115 -45 25],'label','ghi');
% plot_final_dataset('ALOS2_stripe/LOS','fault',fault_file,'lonc',72,'latc',38.5,'ref_lon',71,'misfit_range',30,'max_value',120,'threshold',60,'axis_range',[45 115 -45 25],'label','def');
