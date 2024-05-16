close all
clc
clear

%% setup MATLAB and GMT paths
package_path = '/Users/zej011/work/zeyu/matlab/finite_fault_inversion';
addpath(genpath(package_path));
external_path = '/Users/zej011/work/zeyu/matlab/external_func';
addpath(genpath(external_path));
addpath('/Users/zej011/work/zeyu/matlab');
addpath('/Users/zej011/work/Kang_tutorial/candis');
addpath('/Users/zej011/work/Kang_tutorial/codes_utilities/matlab/igppsar');
% setup your own CSHELL scripts and GMT paths
setenv('PATH',[getenv('PATH'),':/Users/zej011/work/zeyu/cshell']);
setenv('PATH',[getenv('PATH'),':/usr/local/bin']);  % add the path of GMT


%% remove the near-field unwrapping errors manually first
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
% make_insar_data(data_list,Nmin,Nmax,'method','uniform','fault',fault_file,'area',[71.8 73.9 37.7 39.1],'ref_lon',71);
% 
% fault_file = 'all_faults';
% data_list = 'rng.list';         % rng data
% Nmin = 10;    Nmax = 10;
% make_insar_data(data_list,Nmin,Nmax,'method','uniform','fault',fault_file,'area',[71.8 73.9 37.7 39.1],'ref_lon',71);


%% detrending or remove the phase ambiguity
% % build up the coarse model geometry
% W_coarse = 20e3;     % because there is not much slip below 15km depth
% N_layer = 5;
% lp_top = 4e3;
% bias_lp=1.0;
% bias_wp=1.0;
% dip_change_id = [1 2 3];
% % H_offset = [-1.0 -0.2 0.5] * 1e3;   % the segment change the offset from -0.6 to -1.0
% % dip_angle = acosd(H_offset ./ 10e3);
% dip_angle = [87.7, 81.8, 89.3];
% slip_model_vs = load_fault_one_plane('all_faults','width',W_coarse,'layers',N_layer,'len_top',lp_top,'lonc',72,'latc',38.5, ...
%                                      'dip_change_id',dip_change_id,'dip',dip_angle,'l_ratio',bias_lp,'w_ratio',bias_wp,'ref_lon',71);
% slip_model_ds = [];
% 
% % curr_id = max(slip_model_vs(:,1));
% % slip_model_ds = load_fault_dip_shallow('fault_geometry/model_A_splay','len_top',2e3,'fault_id',curr_id, ... 
% %                                        'dip_change_id',[1:7],'dip_depth',[3.2e3,3.2e3,4.2e3,4.2e3,4.2e3,4.2e3,3e3]);
% slip_model = [slip_model_vs;slip_model_ds];
% nflt = max(slip_model(:,1));
% disp(['There are ',num2str(nflt),' fault segments']);
% show_slip_model(slip_model,'ref_lon',71);

% % detrend from the coarse model first
% segment_file = 'seg_connect';
% intersect_file = [];
% detrend_from_inversion(slip_model_vs,slip_model_ds,'detrend.list','ramp_order',[4,4,4], ...
%                       'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file);
% detrend_from_inversion(slip_model_vs,slip_model_ds,'rng.list','ramp_order',[4,4], ...
%                       'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file);

% no GPS data for the Pamir
% just use the benchmark to locate the far-field zero points
% grdin = '/Users/zej011/coseismic/ASC100/LOS/unwrap_clean_sample.grd';
% grdout = '/Users/zej011/coseismic/ASC100/LOS/los_clean_detrend.grd';
% lonf = 73.687830;   latf = 37.860004;  
% ref_lon = 71;   threshold = 0.5;
% remove_ref_from_grid(grdin,grdout,lonf,latf,ref_lon,threshold);

grdin = '/Users/zej011/coseismic/ASC100/LOS2/unwrap_clean_sample.grd';
grdout = '/Users/zej011/coseismic/ASC100/LOS2/los_clean_detrend.grd';
lonf = 74.199446;   latf = 38.064417;  
ref_lon = 71;   threshold = 0.5;
remove_ref_from_grid(grdin,grdout,lonf,latf,ref_lon,threshold);

% grdin = '/Users/zej011/coseismic/DES5/LOS/unwrap_clean_sample.grd';
% grdout = '/Users/zej011/coseismic/DES5/LOS/los_clean_detrend.grd';
% lonf = 73.687830;   latf = 37.860004;  
% ref_lon = 71;   threshold = 0.5;
% remove_ref_from_grid(grdin,grdout,lonf,latf,ref_lon,threshold);

grdin = '/Users/zej011/coseismic/DES5/LOS2/unwrap_clean_sample.grd';
grdout = '/Users/zej011/coseismic/DES5/LOS2/los_clean_detrend.grd';
lonf = 73.215150;   latf = 37.754558;  
ref_lon = 71;   threshold = 0.5;
remove_ref_from_grid(grdin,grdout,lonf,latf,ref_lon,threshold);

% grdin = '/Users/zej011/coseismic/DES5/LOS3/unwrap_clean_sample.grd';
% grdout = '/Users/zej011/coseismic/DES5/LOS3/los_clean_detrend.grd';
% lonf = 73.215150;   latf = 37.754558;  
% ref_lon = 71;   threshold = 0.5;
% remove_ref_from_grid(grdin,grdout,lonf,latf,ref_lon,threshold);

% grdin = '/Users/zej011/coseismic/ALOS2_SCAN/unwrap_clean_sample.grd';
% grdout = '/Users/zej011/coseismic/ALOS2_SCAN/los_clean_detrend.grd';
% lonf = 73.687830;   latf = 37.860004;  
% ref_lon = 71;   threshold = 0.5;
% remove_ref_from_grid(grdin,grdout,lonf,latf,ref_lon,threshold);

% grdin = '/Users/zej011/coseismic/ALOS2_SCAN2/unwrap_clean_sample.grd';
% grdout = '/Users/zej011/coseismic/ALOS2_SCAN2/los_clean_detrend.grd';
% lonf = 73.782137;   latf = 37.679650;  
% ref_lon = 71;   threshold = 0.5;
% remove_ref_from_grid(grdin,grdout,lonf,latf,ref_lon,threshold);

% grdin = '/Users/zej011/coseismic/ALOS2_stripe/LOS/unwrap_clean_sample.grd';
% grdout = '/Users/zej011/coseismic/ALOS2_stripe/LOS/los_clean_detrend.grd';
% lonf = 72.635475;   latf = 38.018589;  
% ref_lon = 71;   threshold = 0.5;
% remove_ref_from_grid(grdin,grdout,lonf,latf,ref_lon,threshold);


%% detrend and clean the offset data                    
% %% detrend the range offsets from LOS data (ASC64 + DES71)
% DEMfile = 'dem_samp.grd';       % all these three files are named uniformly
% losfile = 'los_cut.grd';
% rngfile = 'unwrap_clean_sample.grd';
% curr_dir = pwd;
% % 
% % ASC100_track = '/Users/zej011/coseismic/ASC100/offsets';
% % cd(ASC100_track);
% % !cut_rng_from_los.csh
% % cd(curr_dir);
% % detrend_range_offsets(ASC100_track,losfile,rngfile,DEMfile,'bi_ramp','misfit_range',100,'ref_lon',ref_lon,'lonc',72,'latc',38.5);
% 
% DES5_track = '/Users/zej011/coseismic/DES5/offsets';
% cd(DES5_track);
% !cut_rng_from_los.csh
% cd(curr_dir);
% detrend_range_offsets(DES5_track,losfile,rngfile,DEMfile,'bi_ramp','misfit_range',100,'ref_lon',ref_lon,'lonc',72,'latc',38.5);
% 
% ALOS2_stripe = '/Users/zej011/coseismic/ALOS2_stripe/range';
% cd(ALOS2_stripe);
% !cut_rng_from_los.csh
% cd(curr_dir);
% detrend_range_offsets(ALOS2_stripe,losfile,rngfile,DEMfile,'bi_ramp','misfit_range',100,'ref_lon',ref_lon,'lonc',72,'latc',38.5);

% %% apply the sign mask for the detrended offset data
% fid = fopen('ASC100_test');
% C = textscan(fid,'%s  %d\n');
% data_path = C{:,1};     ndata = length(data_path);
% fclose(fid);
% curr_dir = pwd;
% 
% for jj = 1:ndata
%     this_track = data_path{jj};
%     cd(this_track);
%     movefile los_clean_detrend.grd los_clean_unmask.grd
%     cd(curr_dir);
%     sign_mask_offset(this_track,'los_clean_unmask.grd');
% end
