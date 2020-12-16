function [slip_model,RMS_misfit,model_roughness] = make_fault_from_insar(slip_model_vs,slip_model_ds,iter_step,varargin)
% Build the finite fault model using fault trace derived from both offsets and seismicity data
% return the variance reduction between the model and data
% Started by Zeyu Jin on 07/15/2019

    %% default values
    lambda = 1e-1;
    alpha = 0.25;  % relative weight of RNG data
    beta = 1;   % relative weight of ALOS-2 data
    gamma = 0.8;  % relative weight of AZO data (from CSK)
    segment_file = [];    intersect_file = [];
    shallow_dip_id = [];
    model_type = 'okada';
    fault_file = [];
    ref_lon = 71;
    lon_eq = -117.5;
    lat_eq = 35.5;
    
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
                    case 'ref_lon'
                        ref_lon = varargin{CC*2};
                    case 'lonc'
                        lon_eq = varargin{CC*2};
                    case 'latc'
                        lat_eq = varargin{CC*2};
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
%     [G1_raw,G1,bd1_raw,bd1] = build_green_function(slip_model,['ERS170/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
%     [G2_raw,G2,bd2_raw,bd2] = build_green_function(slip_model,['ERS442/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
%     [G3_raw,G3,bd3_raw,bd3] = build_green_function(slip_model,['ERS120/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);

    [G1_raw,G1,bd1_raw,bd1] = build_green_function(slip_model,['ASC100/LOS/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
    [G2_raw,G2,bd2_raw,bd2] = build_green_function(slip_model,['DES5/LOS/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
%     [G3_raw,G3,bd3_raw,bd3] = build_green_function(slip_model,['ALOS2/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
    G3_raw = []; G3 = []; bd3_raw = []; bd3 = [];
    G4_raw = []; G4 = []; bd4_raw = []; bd4 = [];
    G5_raw = []; G5 = []; bd5_raw = []; bd5 = [];
%     [G4_raw,G4,bd4_raw,bd4] = build_green_function(slip_model,['ALOS2_stripe/LOS/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
%     [G5_raw,G5,bd5_raw,bd5] = build_green_function(slip_model,['ALOS2_stripe/MAI2/los_samp',num2str(iint),'.mat'],'AZO','noramp',model_type);
    
%     [G1_raw,G1,bd1_raw,bd1] = build_green_function(slip_model,['ASC64/branch_cut/three_subswath/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
%     [G2_raw,G2,bd2_raw,bd2] = build_green_function(slip_model,['DES71/branch_cut/three_subswath/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
%     [G3_raw,G3,bd3_raw,bd3] = build_green_function(slip_model,['ALOS-2/T065/three_subswath/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
%     [G4_raw,G4,bd4_raw,bd4] = build_green_function(slip_model,['ALOS-2/T066/three_subswath/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
%     [G5_raw,G5,bd5_raw,bd5] = build_green_function(slip_model,['ASC64/offsets/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
%     [G6_raw,G6,bd6_raw,bd6] = build_green_function(slip_model,['DES71/offsets/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
%     [G7_raw,G7,bd7_raw,bd7] = build_green_function(slip_model,['Cosmo_Skymed/ASC_offsets/los_samp',num2str(iint),'.mat'],'AZO','noramp',model_type);
%     [G8_raw,G8,bd8_raw,bd8] = build_green_function(slip_model,['Cosmo_Skymed/DES_offsets/los_samp',num2str(iint),'.mat'],'AZO','noramp',model_type);
% %     [Gp_raw,Gp,bp_raw,bp] = build_green_function(slip_model,'GPS/continue_GPS/continuous_gps_3d.mat','cgps','noramp',model_type,0.3);     % default weight of cGPS is 0.3
% %     [Gg_raw,Gg,bg_raw,bg] = build_green_function(slip_model,'GPS/continue_GPS/Peng_far_GPS/Peng_far_gps_3d.mat','cgps','noramp',model_type,0.3);
% %     [Gp_raw,Gp,bp_raw,bp] = build_green_function(slip_model,'GPS/Floyd_GPS/gps_syn.mat','camp_gps','noramp',model_type,0.3);
% %     [Gs_raw,Gs,bs_raw,bs] = build_green_function(slip_model,'GPS/survey_GPS/gps_syn.mat','camp_gps','noramp',model_type,0.25);
%     [Gp_raw,Gp,bp_raw,bp] = build_green_function(slip_model,'GPS/Floyd_GPS/Floyd_GPS_all.mat','camp_gps','noramp',model_type,0.3);
%     [Gs_raw,Gs,bs_raw,bs] = build_green_function(slip_model,'GPS/survey_GPS/survey_gps_2d.mat','camp_gps','noramp',model_type,0.25);


    %% generate Green's function and smooth matrix
    [H,h1,~] = build_smooth_function(slip_model_vs,slip_model_ds,segment_file,intersect_file,'noramp','dip_id',shallow_dip_id); 
    
% %     segment_ID = 1:2;  top_layer_no = 1; ratio = 2e-4;
%     segment_ID = 3;  top_layer_no = 1; ratio = 2e-4;
%     [Ws,ds] = zero_slip_boundary(slip_model,segment_ID,top_layer_no,ratio);    % for the top layer
%     
%     plane_fault = 1:5; 
    plane_fault = 1:4;
%     plane_fault = 1:2;
    bottom_layer_no = max(slip_model(:,3)); ratio = 5e-4;
    [Wb,db] = zero_slip_boundary(slip_model,plane_fault,bottom_layer_no,ratio);    % for the bottom layer
%       Wb = []; db = [];
% %     plane_fault = 2; ratio = 6e-4;
% %     [Wb2,db2] = zero_slip_boundary(slip_model,plane_fault,bottom_layer_no,ratio);    % for the bottom layer 
% %     Wb = [Wb1;Wb2];
% %     db = [db1;db2];
%     plane_fault = 1:8;  bottom_layer_no = max(slip_model(:,3)); ratio = 1e-4;
%     [Wb,db] = zero_slip_boundary(slip_model,plane_fault,bottom_layer_no,ratio);    % for the bottom layer
%     
% %     left_fault = 1;  ratio = 3e-4;
% %     [Wl1,dl1] = zero_slip_boundary(slip_model,left_fault,'left',ratio);   
% %     left_fault = 2;  ratio = 6e-4;
% %     [Wl2,dl2] = zero_slip_boundary(slip_model,left_fault,'left',ratio);
% %     Wl = [Wl1;Wl2];
% %     dl = [dl1;dl2]; 
%     left_fault = 4:5;  
    left_fault = 3:4;
%     left_fault = 2;
    ratio = 3e-4;
    [Wl,dl] = zero_slip_boundary(slip_model,left_fault,'left',ratio);
%       Wl = []; dl = [];
%     
% %     right_fault = 1;   ratio = 5e-4;
% %     [Wr1,dr1] = zero_slip_boundary(slip_model,right_fault,'right',ratio);
% %     right_fault = 2;   ratio = 6e-4;
% %     [Wr2,dr2] = zero_slip_boundary(slip_model,right_fault,'right',ratio);   
% %     Wr = [Wr1;Wr2];
% %     dr = [dr1;dr2];
    right_fault = 1;   ratio = 3e-4;
    [Wr,dr] = zero_slip_boundary(slip_model,right_fault,'right',ratio);
%       Wr = [];  dr = [];
%     
% %     % constrained by moment (M5.8)
% %     fault_id = slip_model(:,1);
% %     indx_M54 = fault_id == 1;
% %     indx_M58 = fault_id == 2;
% %     Mw2 = 5.8;
% %     [area,M0] = constrain_by_moment(slip_model,model_type,Mw2,indx_M58);
% %     RR = 1.8e-4;
% % %     area = [];
% % %     M0 = [];

    % adjust the relative weight between ASC and DES track to be 1:1
    % better for fitting one track and decomposition
%     A2D = 1/3;    % three ascending tracks (ASC64,T065,T066) & one descending track (DES71)
%     A2D = 1;       % equal weight of each dataset
%     G_raw = [A2D*G1;G2;A2D*beta*G3;A2D*beta*G4;alpha*G5;alpha*G6;gamma*G7;gamma*G8;Gp;Gs];  
%     G_rough = [H*lambda/h1;Ws;Wb;Wl;Wr];
% %     G_rough = [Ws;Wb;Wl;Wr];
    G_raw = [G1;beta*G2;gamma*G3;G4;G5]; 
    Greens = [G_raw;H*lambda/h1;Wb;Wl;Wr];
    
%     bd_raw = [A2D*bd1;bd2;A2D*beta*bd3;A2D*beta*bd4;alpha*bd5;alpha*bd6;gamma*bd7;gamma*bd8;bp;bs];
%     bdata_sm = [bd_raw;zeros(h1,1);ds;db;dl;dr];
% %     bdata_sm = [bd_raw;ds;db;dl;dr];
    bd_raw = [bd1;beta*bd2;gamma*bd3;bd4;bd5];
    bdata_sm = [bd_raw;zeros(h1,1);db;dl;dr];
    
%     % for 1995 M5+ events
%     G_raw = [G1;G2;G3];
%     Greens = [G_raw;H*lambda/h1;Ws;Wb;Wl;Wr;RR*area];
%     bd_raw = [bd1;bd2;bd3];
%     bdata_sm = [bd_raw;zeros(h1,1);ds;db;dl;dr;RR*M0];
    
%     GrF = [G1_raw;G2_raw;G3_raw;G4_raw;G5_raw;G6_raw;G7_raw;G8_raw;Gp_raw;Gs_raw];   % without any weight
%     Bdata = [bd1_raw;bd2_raw;bd3_raw;bd4_raw;bd5_raw;bd6_raw;bd7_raw;bd8_raw;bp_raw;bs_raw];
% %     GrF = [G1_raw;G2_raw;G3_raw];
% %     Bdata = [bd1_raw;bd2_raw;bd3_raw];
    GrF = [G1_raw;G2_raw;G3_raw;G4_raw;G5_raw];
    Bdata = [bd1_raw;bd2_raw;bd3_raw;bd4_raw;bd5_raw];

    %% the postivity constraint (adapted from Yuri's code)
    nflt = max(slip_model(:,1));
    tSm = zeros(1,nflt+1);
    fault_id = slip_model(:,1);
    for i=1:nflt
        tSm(i+1) = length(find(fault_id == i));
    end
    add_col = 0; NT = 2; NS = nflt;  % the number of segments  
%     [lb,ub] = bounds_new_C(NS,NT,tSm,add_col,slip_model);
%     [lb,ub] = bounds_new_B(NS,NT,tSm,add_col,slip_model);
%     [lb,ub] = bounds_new_A(NS,NT,tSm,add_col);
    [lb,ub] = bounds_new(NS,NT,tSm,add_col);
%     [lb,ub] = bounds_new_M5(NS,NT,tSm,add_col);
%     [lb,ub] = bounds_resolution(NS,NT,tSm,add_col);

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
    slip_model(:,12) = u(1:sum(tSm));
    slip_model(:,13) = u(sum(tSm)+1:end);
%     show_slip_model(slip_model,'seismicity_profile/Ross_seismicity_cut.mat');
%     show_slip_model(slip_model,'misfit_range',400,'ref_lon',ref_lon,'lonc',lon_eq,'latc',lat_eq,'axis_range',[60 110 -45 25 -25 0]);
    show_slip_model(slip_model,'misfit_range',25,'ref_lon',ref_lon,'lonc',lon_eq,'latc',lat_eq,'axis_range',[80 120 -10 40 -30 0]);
%     show_slip_model(slip_model,'misfit_range',25,'ref_lon',ref_lon,'lonc',lon_eq,'latc',lat_eq,'axis_range',[80 120 -10 40 -12 0]);
%     show_slip_model(slip_model,'axis_range',[-25 0 15 45 -15 0],'misfit_range',10);
    write_slip_model_okada(slip_model,'fault_M7.slip');
%     write_slip_model_okada(slip_model,'fault_M5.slip');

    %% plot the resampled data fitting
    insar_model1 = G1_raw * u;
    insar_model2 = G2_raw * u;
%     insar_model3 = G3_raw * u;
%     insar_model4 = G4_raw * u;
%     insar_model5 = G5_raw * u;
%     insar_model1 = G1_raw * u;
%     insar_model2 = G2_raw * u;
%     insar_model3 = G3_raw * u;
%     insar_model4 = G4_raw * u;
%     insar_model5 = G5_raw * u;
%     insar_model6 = G6_raw * u;
%     insar_model7 = G7_raw * u;
%     insar_model8 = G8_raw * u;
% %     cgps_model = Gg_raw * u; 
%     Floyd_gps_model = Gp_raw * u;
%     survey_gps_model = Gs_raw * u;  
    
    
    plot_insar_model_resampled(['ASC100/LOS/los_samp',num2str(iint),'.mat'],insar_model1,'iter_step',iint,'fault',fault_file,'model_type',model_type,'misfit_range',8,'defo_max',8,'ref_lon',ref_lon,'lonc',lon_eq,'latc',lat_eq);
    plot_insar_model_resampled(['DES5/LOS/los_samp',num2str(iint),'.mat'],insar_model2,'iter_step',iint,'fault',fault_file,'model_type',model_type,'misfit_range',8,'defo_max',8,'ref_lon',ref_lon,'lonc',lon_eq,'latc',lat_eq);
%     plot_insar_model_resampled(['ALOS2/los_samp',num2str(iint),'.mat'],insar_model3,'iter_step',iint,'fault',fault_file,'model_type',model_type,'misfit_range',8,'defo_max',8,'ref_lon',ref_lon,'lonc',lon_eq,'latc',lat_eq);
%     plot_insar_model_resampled(['ALOS2_stripe/LOS/los_samp',num2str(iint),'.mat'],insar_model4,'iter_step',iint,'fault',fault_file,'model_type',model_type,'misfit_range',30,'ref_lon',ref_lon,'lonc',lon_eq,'latc',lat_eq);
%     plot_insar_model_resampled(['ALOS2_stripe/MAI2/los_samp',num2str(iint),'.mat'],insar_model5,'iter_step',iint,'fault',fault_file,'model_type',model_type,'misfit_range',40,'ref_lon',ref_lon,'lonc',lon_eq,'latc',lat_eq);
    
%     plot_insar_model_resampled(['ASC64/branch_cut/three_subswath/los_samp',num2str(iint),'.mat'],insar_model1,'iter_step',iint,'fault',fault_file,'axis_range',[-60,60,-20,70],'model_type',model_type);
%     plot_insar_model_resampled(['DES71/branch_cut/three_subswath/los_samp',num2str(iint),'.mat'],insar_model2,'iter_step',iint,'fault',fault_file,'axis_range',[-55,50,-20,75],'model_type',model_type);
%     plot_insar_model_resampled(['ALOS-2/T065/three_subswath/los_samp',num2str(iint),'.mat'],insar_model3,'iter_step',iint,'fault',fault_file,'axis_range',[-35,35,-5,60],'model_type',model_type);
%     plot_insar_model_resampled(['ALOS-2/T066/three_subswath/los_samp',num2str(iint),'.mat'],insar_model4,'iter_step',iint,'fault',fault_file,'axis_range',[-45,35,-10,60],'model_type',model_type);
%     plot_insar_model_resampled(['ASC64/offsets/los_samp',num2str(iint),'.mat'],insar_model5,'misfit_range',20,'iter_step',iint,'fault',fault_file,'axis_range',[-30,20,0,50],'model_type',model_type);
%     plot_insar_model_resampled(['DES71/offsets/los_samp',num2str(iint),'.mat'],insar_model6,'misfit_range',20,'iter_step',iint,'fault',fault_file,'axis_range',[-30,20,0,50],'model_type',model_type);
%     plot_insar_model_resampled(['Cosmo_Skymed/ASC_offsets/los_samp',num2str(iint),'.mat'],insar_model7,'misfit_range',20,'iter_step',iint,'fault',fault_file,'axis_range',[-30,20,0,50],'model_type',model_type);
%     plot_insar_model_resampled(['Cosmo_Skymed/DES_offsets/los_samp',num2str(iint),'.mat'],insar_model8,'misfit_range',20,'iter_step',iint,'fault',fault_file,'axis_range',[-30,20,0,50],'model_type',model_type);
% %     plot_gps_model('GPS/Floyd_GPS/gps_syn.mat',Floyd_gps_model,'iter_step',iint,'data_type','survey','model_type',model_type);
% %     plot_gps_model('GPS/survey_GPS/gps_syn.mat',survey_gps_model,'iter_step',iint,'data_type','survey','model_type',model_type);
% %     plot_gps_model('GPS/continue_GPS/Peng_far_GPS/Peng_far_gps_3d.mat',cgps_model,'iter_step',iint,'data_type','cont','model_type',model_type);
%    plot_gps_model('GPS/Floyd_GPS/Floyd_GPS_all.mat',Floyd_gps_model,'iter_step',iint,'data_type','survey','model_type',model_type);
%    plot_gps_model('GPS/survey_GPS/survey_gps_2d.mat',survey_gps_model,'iter_step',iint,'data_type','survey','site_name','GPS/survey_GPS/campaign_sites','model_type',model_type);
 
%     % compute the moment contributed by M6 and M7 events individually    
%     all_indx = [1:size(slip_model,1)]';
%     segID = slip_model(:,1);
%     LL = [7 8];       % assume M6.4 slip concentrated on two LL faults
%     M6_indx = find(segID == LL(1) | segID == LL(2));
%     M7_indx = setdiff(all_indx,M6_indx);
%     M6_model = slip_model(M6_indx,:);
%     M7_model = slip_model(M7_indx,:);
%     compute_moment(M6_model,model_type);
%     compute_moment(M7_model,model_type);
%     compute_moment(slip_model,model_type);
    
%     model54 = slip_model(indx_M54,:); 
%     model58 = slip_model(indx_M58,:);
%     compute_moment(model54,model_type);
%     compute_moment(model58,model_type);
      compute_moment(slip_model,model_type);
    
end  