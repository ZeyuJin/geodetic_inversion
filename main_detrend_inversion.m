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

% % remove the near-field unwrapping errors manually first
% clean_insar_data;

% % cut the grid of looking angles and DEM to be the same size
% fid = fopen('all_data_list');
% C = textscan(fid,'%s\n');
% data_path = C{:,1};     ndata = length(data_path);
% fclose(fid);
% curr_dir = pwd;
% for ii = 1:ndata
%     this_path = data_path{ii};
%     cd(this_path);
%     !subsample_grd.csh
%     cd(curr_dir);
% end

% % uniform sampling LOS and AZO data (not apply for RNG yet)
% data_list = 'detrend.list';     % LOS data
% fault_file = 'all_faults';
% Nmin = 20;    Nmax = 20;
% make_insar_data(data_list,Nmin,Nmax,'method','uniform','fault',fault_file);
% 
% fault_file = 'all_faults';
% data_list = 'CSK.list';         % AZO data
% Nmin = 5;    Nmax = 5;
% make_insar_data(data_list,Nmin,Nmax,'method','uniform','fault',fault_file);


% % build up the coarse model geometry
% W_coarse = 15e3;     % because there is not much slip below 15km depth
% N_layer = 5;
% lp_top = 4e3;
% bias_lp=1.0;
% bias_wp=1.0;
% dip_change_id = [1 2 7 8 9 10];
% H_offset = [-0.6 -0.2 0.5 0.5 0.6 -0.4] * 1e3;
% dip_angle = acosd(H_offset ./ 10e3);
% slip_model_vs = load_fault_one_plane('fault_geometry/model_B_plane','width',W_coarse,'layers',N_layer,'len_top',lp_top, ...
%                                      'dip_change_id',dip_change_id,'dip',dip_angle,'l_ratio',bias_lp,'w_ratio',bias_wp);
% % added for model B                                 
% ID = 7;  N_layer_cut = 1;
% [~,nL] = smoo1_each_plane(slip_model_vs);
% patch_this_ID = nL(ID,:);
% total_patch_before = sum(sum(nL(1:ID-1,:),2));
% patch_3_layer_this_ID = 1:sum(patch_this_ID(1:N_layer_cut));
% patch_disappear = total_patch_before + patch_3_layer_this_ID;
% patch_rest_layer_this_ID = max(patch_disappear)+1:sum(sum(nL(1:ID,:),2));
% % recompute the patch number and layer
% slip_model_vs(patch_rest_layer_this_ID,2) = slip_model_vs(patch_rest_layer_this_ID,2) - length(patch_disappear);
% slip_model_vs(patch_rest_layer_this_ID,3) = slip_model_vs(patch_rest_layer_this_ID,3) - N_layer_cut;
% slip_model_vs(patch_disappear,:) = [];                                 
% 
% curr_id = max(slip_model_vs(:,1));
% slip_model_ds = load_fault_dip_shallow('fault_geometry/model_B_splay','len_top',2e3,'fault_id',curr_id,'dip_change_id',[1:4],'dip_depth',[3.2e3,3.2e3,4.2e3,3e3]);
% % slip_model = [slip_model_vs;slip_model_ds];
% % nflt = max(slip_model(:,1));
% % disp(['There are ',num2str(nflt),' fault segments']);
% % show_slip_model(slip_model);
% 
% % detrend from the coarse model first
% segment_file = 'fault_geometry/model_B_seg_coarse';     intersect_file = 'fault_geometry/model_B_dip_coarse';
% detrend_from_inversion(slip_model_vs,slip_model_ds,'detrend.list','ramp_order',[5,4,5,7], ...
%                       'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file);
% detrend_from_inversion(slip_model_vs,slip_model_ds,'CSK.list','data_type','AZO','ramp_order',[4,4], ...
%                       'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file);

                    
% %% detrend the range offsets from LOS data (ASC64 + DES71)
% DEMfile = 'dem_samp.grd';       % all these three files are named uniformly
% losfile = 'los_cut.grd';
% rngfile = 'unwrap_clean_sample.grd';
% curr_dir = pwd;
% 
% ASC64_track = '/Users/zej011/Ridgecrest/data_resample/ASC64/offsets';
% cd(ASC64_track);
% !cut_rng_from_los.csh
% cd(curr_dir);
% detrend_range_offsets(ASC64_track,losfile,rngfile,DEMfile,'sine','misfit_range',30);
% 
% DES71_track = '/Users/zej011/Ridgecrest/data_resample/DES71/offsets';
% cd(DES71_track);
% !cut_rng_from_los.csh
% cd(curr_dir);
% detrend_range_offsets(DES71_track,losfile,rngfile,DEMfile,'qu_ramp_7','misfit_range',20);

% %% apply the sign mask for the detrended offset data
% fid = fopen('offset.list');
% C = textscan(fid,'%s\n');
% data_path = C{:,1};     ndata = length(data_path);
% fclose(fid);
% curr_dir = pwd;
% for jj = 1:ndata
%     this_track = data_path{jj};
%     cd(this_track);
%     movefile los_clean_detrend.grd los_clean_unmask.grd
%     cd(curr_dir);
%     sign_mask_offset(this_track,'los_clean_unmask.grd');
% end
