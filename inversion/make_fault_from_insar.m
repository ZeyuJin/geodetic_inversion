function [slip_model,RMS_misfit,model_roughness] = make_fault_from_insar(slip_model_vs,slip_model_ds,iter_step,varargin)
% Build the finite fault model using fault trace derived from both offsets and seismicity data
% return the variance reduction between the model and data
% Started by Zeyu Jin on 07/15/2019

    %% default values
    lambda = 1e-1;
    alpha = 0.25;  % relative weight of RNG data
    beta = 0.8;   % relative weight of ALOS-2 data
    gamma = 0.1;  % relative weight of AZO data (from CSK)
    segment_file = [];    intersect_file = [];
    shallow_dip_id = [];
    model_type = 'okada';
    fault_file = [];
    
    %% read varargin values and assembly
    if ~isempty(varargin)
        for CC = 1:floor(length(varargin)/2)
            try
                switch lower(varargin{CC*2-1})
                    case 'smoothness'
                        lambda = varargin{CC*2};
                    case 'rng_ratio'
                        alpha = varargin{CC*2};
                    case 'alos_ratio'
                        beta = varargin{CC*2};
                    case 'azo_ratio'
                        gamma = varargin{CC*2};
                    case 'segment_smooth_file'
                        segment_file = varargin{CC*2};
                    case 'intersect_smooth_file'
                        intersect_file = varargin{CC*2};
                    case 'shallow_dip_id'
                        shallow_dip_id = varargin{CC*2};  % to control the dip slip component
                    case 'model_type'
                        model_type = varargin{CC*2};  % homogenous or layered
                    case 'fault'
                        fault_file = varargin{CC*2};    
                end
            catch
                error('Unrecognized Keyword');
            end
        end
    end    

    %% read downsampled data
    iint = iter_step;  % the number of resampling iteration
    slip_model = [slip_model_vs;slip_model_ds];   
    slip_model(:,2)=[1:size(slip_model,1)]';    % recomputed finally to combine all the fault segments
    disp(['There are total ',num2str(max(slip_model(:,1))),' segments']);
    
    [G1_raw,G1,bd1_raw,bd1] = build_green_function(slip_model,['ASC64/branch_cut/three_subswath/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
    [G2_raw,G2,bd2_raw,bd2] = build_green_function(slip_model,['DES71/branch_cut/three_subswath/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
    [G3_raw,G3,bd3_raw,bd3] = build_green_function(slip_model,['ALOS-2/T065/three_subswath/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
    [G4_raw,G4,bd4_raw,bd4] = build_green_function(slip_model,['ALOS-2/T066/three_subswath/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
    [G5_raw,G5,bd5_raw,bd5] = build_green_function(slip_model,['ASC64/offsets/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
    [G6_raw,G6,bd6_raw,bd6] = build_green_function(slip_model,['DES71/offsets/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
    [G7_raw,G7,bd7_raw,bd7] = build_green_function(slip_model,['Cosmo_Skymed/ASC_offsets/los_samp',num2str(iint),'.mat'],'AZO','noramp',model_type);
    [G8_raw,G8,bd8_raw,bd8] = build_green_function(slip_model,['Cosmo_Skymed/DES_offsets/los_samp',num2str(iint),'.mat'],'AZO','noramp',model_type);
    [Gp_raw,Gp,bp_raw,bp] = build_green_function(slip_model,'GPS/continue_GPS/continuous_gps_3d.mat','cgps','noramp',model_type,0.3);     % default weight of cGPS is 0.3
    [Gs_raw,Gs,bs_raw,bs] = build_green_function(slip_model,'GPS/survey_GPS/survey_gps_2d.mat','camp_gps','noramp',model_type,0.25);
%     Gs = [];  bs = [];
%     Gs_raw = [];  bs_raw = [];   % test using campaign GPS data or not

    %% generate Green's function and smooth matrix
    [H,h1,~] = build_smooth_function(slip_model_vs,slip_model_ds,segment_file,intersect_file,'noramp','dip_id',shallow_dip_id); 
    segment_ID = 3;  top_layer_no = 1; ratio = 0.01;
    [Ws,ds] = zero_slip_boundary(slip_model,segment_ID,top_layer_no,ratio);
%     Ws = [];   ds = [];
    % adjust the relative weight between ASC and DES track to be 1:1
    % better for fitting one track and decomposition
%     A2D = 1/3;    % three ascending tracks (ASC64,T065,T066) & one descending track (DES71)
    A2D = 1;       % equal weight of each dataset
    
    G_raw = [A2D*G1;G2;A2D*beta*G3;A2D*beta*G4;alpha*G5;alpha*G6;gamma*G7;gamma*G8;Gp;Gs];  
    Greens = [G_raw;H*lambda/h1;Ws];
    bd_raw = [A2D*bd1;bd2;A2D*beta*bd3;A2D*beta*bd4;alpha*bd5;alpha*bd6;gamma*bd7;gamma*bd8;bp;bs]; 
    bdata_sm = [bd_raw;zeros(h1,1);ds];
    
    GrF = [G1_raw;G2_raw;G3_raw;G4_raw;G5_raw;G6_raw;G7_raw;G8_raw;Gp_raw;Gs_raw];   % without any weight
    Bdata = [bd1_raw;bd2_raw;bd3_raw;bd4_raw;bd5_raw;bd6_raw;bd7_raw;bd8_raw;bp_raw;bs_raw];

    %% the postivity constraint (adapted from Yuri's code)
    nflt = max(slip_model(:,1));
    tSm = zeros(1,nflt+1);
    fault_id = slip_model(:,1);
    for i=1:nflt
        tSm(i+1) = length(find(fault_id == i));
    end
    add_col = 0; NT = 2; NS = nflt;  % the number of segments  
%    [lb,ub] = bounds_new_A(NS,NT,tSm,add_col,slip_model);
    [lb,ub] = bounds_new_A(NS,NT,tSm,add_col,slip_model);

    % linear inversion
    options = optimset('LargeScale','on','DiffMaxChange',1e-1,'DiffMinChange',1e-12, ...
        'TolCon',1e-12,'TolFun',1e-12,'TolPCG',1e-12,'TolX',1e-12,'MaxIter',1e9,'MaxPCGIter',1e9);
    [u,resnorm,residual,exitflag] = lsqlin(Greens,double(bdata_sm),[],[],[],[],lb,ub,[],options);
    
    % compute the reduction of total variance (before weighing) of the downsampled data
    rms0 = sum(Bdata.^2);
    rms = sum((GrF*u-Bdata).^2);
    redu_perc = 100*(rms0-rms)/rms0;   
    fprintf('rms misfit (dat., res.) = %e %e (%f%%) \n',rms0,rms,redu_perc);
    fprintf('resnorm, resid. = %e %e \n',sqrt(resnorm),mean(residual));
    fprintf('exitflag is %d\n',exitflag);   % 1 means the function converged to a solution x
    
    rough_matrix = H*u; 
    RMS_misfit = sum((G_raw*u - bd_raw).^2);   % use chi-square statistic instead
    model_roughness = sqrt(sum(rough_matrix.^2)/length(rough_matrix));


    %% plot and save the finite fault inversion
%     slip_model(:,12) = u(1:sum(tSm));
%     slip_model(:,13) = u(sum(tSm)+1:end);
% %     show_slip_model(slip_model,'seismicity_profile/Ross_seismicity_cut.mat');
%     show_slip_model(slip_model,'misfit_range',500);
% %     write_slip_model_okada(slip_model,'fault_M7.slip');
% 
%     %% plot the resampled data fitting
%     insar_model1 = G1_raw * u;
%     insar_model2 = G2_raw * u;
%     insar_model3 = G3_raw * u;
%     insar_model4 = G4_raw * u;
%     insar_model5 = G5_raw * u;
%     insar_model6 = G6_raw * u;
%     insar_model7 = G7_raw * u;
%     insar_model8 = G8_raw * u;
%     cgps_model = Gp_raw * u;     
%     survey_gps_model = Gs_raw * u;    
%     
%     plot_insar_model_resampled(['ASC64/branch_cut/three_subswath/los_samp',num2str(iint),'.mat'],insar_model1,'iter_step',iint,'fault',fault_file,'axis_range',[-60,60,-20,70],'model_type',model_type);
%     plot_insar_model_resampled(['DES71/branch_cut/three_subswath/los_samp',num2str(iint),'.mat'],insar_model2,'iter_step',iint,'fault',fault_file,'axis_range',[-55,50,-20,75],'model_type',model_type);
%     plot_insar_model_resampled(['ALOS-2/T065/three_subswath/los_samp',num2str(iint),'.mat'],insar_model3,'iter_step',iint,'fault',fault_file,'axis_range',[-35,35,-5,60],'model_type',model_type);
%     plot_insar_model_resampled(['ALOS-2/T066/three_subswath/los_samp',num2str(iint),'.mat'],insar_model4,'iter_step',iint,'fault',fault_file,'axis_range',[-45,35,-10,60],'model_type',model_type);
%     plot_insar_model_resampled(['ASC64/offsets/los_samp',num2str(iint),'.mat'],insar_model5,'misfit_range',50,'iter_step',iint,'fault',fault_file,'axis_range',[-30,20,0,50],'model_type',model_type);
%     plot_insar_model_resampled(['DES71/offsets/los_samp',num2str(iint),'.mat'],insar_model6,'misfit_range',50,'iter_step',iint,'fault',fault_file,'axis_range',[-30,20,0,50],'model_type',model_type);
%     plot_insar_model_resampled(['Cosmo_Skymed/ASC_offsets/los_samp',num2str(iint),'.mat'],insar_model7,'misfit_range',50,'iter_step',iint,'fault',fault_file,'axis_range',[-30,20,0,50],'model_type',model_type);
%     plot_insar_model_resampled(['Cosmo_Skymed/DES_offsets/los_samp',num2str(iint),'.mat'],insar_model8,'misfit_range',50,'iter_step',iint,'fault',fault_file,'axis_range',[-30,20,0,50],'model_type',model_type);
%     plot_gps_model('GPS/continue_GPS/continuous_gps_3d.mat',cgps_model,'iter_step',iint,'data_type','cont','model_type',model_type);
%     plot_gps_model('GPS/survey_GPS/survey_gps_2d.mat',survey_gps_model,'iter_step',iint,'data_type','survey','site_name','GPS/survey_GPS/campaign_sites','model_type',model_type);
%     
%     % compute the scalar seismic moment
%     strike_u = slip_model(:,12) ./ 100;     % in meters
%     strike_d = slip_model(:,13) ./ 100;
%     dz_top = slip_model(:,6) ./ 1000;
%     dip = slip_model(:,10);
%     D = sqrt(strike_u.^2 + strike_d.^2);
%     lpatch = slip_model(:,7);
%     wpatch = slip_model(:,8);
%     Apatch = lpatch .* wpatch;   
%     patch_center = dz_top*(-1) + wpatch.*sind(dip)/2/1000;
%     
%     depth = [0 3 6 9 12 15 20 29];
%     shear = [2.32 2.71 3.06 3.21 3.35 3.52 3.76 4.25] * 10e9; 
%     
%     if strcmp(model_type,'okada')
%         mu = 30e9;
%         M0 = sum(mu .* D .* Apatch);
%         Mw = 2/3*(log10(M0) - 9.1);
%         disp(['The moment magnitude is Mw = ',num2str(Mw)]);
%         fprintf('\n');
%     else
%         total = 0;
%         for kk = 1:length(patch_center)
%             for mm = 2:length(depth)
%                 if depth(mm) > patch_center(kk)
%                     this_depth = depth(mm);
%                     top_depth = depth(mm-1);
%                     bottom_depth = depth(mm);
%                     top_shear = shear(mm-1);
%                     bottom_shear = shear(mm);
%                     mu_this_depth = (top_shear*(bottom_depth-this_depth)+bottom_shear*(this_depth-top_depth))/(bottom_depth-top_depth);
%                     tmp = D(kk) * Apatch(kk) * mu_this_depth;  
%                     break
%                 end
%             end
%         total = total + tmp;
%         end
%         Mw = 2/3*(log10(total) - 9.1);
%         disp(['The moment magnitude is Mw = ',num2str(Mw)]);
%         fprintf('\n');
%     end

end  