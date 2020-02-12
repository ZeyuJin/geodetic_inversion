function dataset_sensitivity(test_data_type,slip_model_vs,slip_model_ds,iter_step,varargin)
    %% default values
    lambda = 1e-1;
    segment_file = [];    intersect_file = [];
    shallow_dip_id = [];
    model_type = 'okada';
    alpha = [4:13] * 0.1;  % for ALOS-2 and cGPS
%     alpha = [1:10] * 0.05;   % for S1A-RNG and CSK-AZO
%     alpha = [2:11] * 0.1;    % for survey GPS
    
    %% read varargin values and assembly
    if ~isempty(varargin)
        for CC = 1:floor(length(varargin)/2)
            try
                switch lower(varargin{CC*2-1})
                    case 'smoothness'
                        lambda = varargin{CC*2};
                    case 'segment_smooth_file'
                        segment_file = varargin{CC*2};
                    case 'intersect_smooth_file'
                        intersect_file = varargin{CC*2};
                    case 'shallow_dip_id'
                        shallow_dip_id = varargin{CC*2};  % to control the dip slip component
                    case 'model_type'
                        model_type = varargin{CC*2};  % homogenous or layered
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
    
    % use S1A InSAR data as the reference
    % all the other data use relative weighting to S1A InSAR data
    [G1_raw,G1,bd1_raw,bd1] = build_green_function(slip_model,['ASC64/branch_cut/three_subswath/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
    [G2_raw,G2,bd2_raw,bd2] = build_green_function(slip_model,['DES71/branch_cut/three_subswath/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
    S1A_GrF = [G1;G2];
    S1A_Bdata = [bd1;bd2];
    
    if strcmp(test_data_type,'ALOS2')
        [G3_raw,G3,bd3_raw,bd3] = build_green_function(slip_model,['ALOS-2/T065/three_subswath/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
        [G4_raw,G4,bd4_raw,bd4] = build_green_function(slip_model,['ALOS-2/T066/three_subswath/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
    elseif strcmp(test_data_type,'RNG')
        [G3_raw,G3,bd3_raw,bd3] = build_green_function(slip_model,['ASC64/offsets/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
        [G4_raw,G4,bd4_raw,bd4] = build_green_function(slip_model,['DES71/offsets/los_samp',num2str(iint),'.mat'],'insar','noramp',model_type);
    elseif strcmp(test_data_type,'AZO')
        [G3_raw,G3,bd3_raw,bd3] = build_green_function(slip_model,['Cosmo_Skymed/ASC_offsets/los_samp',num2str(iint),'.mat'],'AZO','noramp',model_type);
        [G4_raw,G4,bd4_raw,bd4] = build_green_function(slip_model,['Cosmo_Skymed/DES_offsets/los_samp',num2str(iint),'.mat'],'AZO','noramp',model_type);
    elseif strcmp(test_data_type,'cgps')
        [G3_raw,G3,bd3_raw,bd3] = build_green_function(slip_model,'GPS/continue_GPS/continuous_gps_3d.mat','cgps','noramp',model_type,1); 
        G4_raw = [];  G4 = [];  bd4_raw = [];  bd4 = []; 
    elseif strcmp(test_data_type,'camp_gps')
        [G3_raw,G3,bd3_raw,bd3] = build_green_function(slip_model,'GPS/survey_GPS/survey_gps_2d.mat','camp_gps','noramp',model_type,1);
        G4_raw = [];  G4 = [];  bd4_raw = [];  bd4 = [];
    else
        error('There is something wrong with the input data type!');
    end
    test_GrF = [G3;G4];
    test_Bdata = [bd3;bd4];
    
%     GrF = [G1_raw;G2_raw;G3_raw;G4_raw];
%     Bdata = [bd1_raw;bd2_raw;bd3_raw;bd4_raw];

    %% generate Green's function and smooth matrix
    [H,h1,~] = build_smooth_function(slip_model_vs,slip_model_ds,segment_file,intersect_file,'noramp','dip_id',shallow_dip_id); 
    
    %% the postivity constraint (adapted from Yuri's code)
    nflt = max(slip_model(:,1));
    tSm = zeros(1,nflt+1);
    fault_id = slip_model(:,1);
    for i=1:nflt
        tSm(i+1) = length(find(fault_id == i));
    end
    add_col = 0; NT = 2; NS = nflt;  % the number of segments  
%     fault_patch = slip_model(:,1:3);
%     [lb,ub] = bounds_new(NS,NT,tSm,add_col,fault_patch,smooth_patch);  
    [lb,ub] = bounds_new_A(NS,NT,tSm,add_col,slip_model);
    
    %% loop the relative weighting
    L_curve = zeros(length(alpha),2);
    for ii = 1:length(alpha)
        ratio = alpha(ii);
        G_raw = [G1;G2;ratio*G3;ratio*G4];        Greens = [G_raw;H*lambda/h1];
        bd_raw = [bd1;bd2;ratio*bd3;ratio*bd4];   bdata_sm = [bd_raw;zeros(h1,1)];
        
        % linear inversion
        options = optimset('LargeScale','on','DiffMaxChange',1e-1,'DiffMinChange',1e-12, ...
            'TolCon',1e-12,'TolFun',1e-12,'TolPCG',1e-12,'TolX',1e-12,'MaxIter',1e9,'MaxPCGIter',1e9);
        [u,resnorm,residual,exitflag] = lsqlin(Greens,double(bdata_sm),[],[],[],[],lb,ub,[],options); 
       
        % compute the reduction of total variance (before weighing) of the downsampled data
%         rms0 = sum(Bdata.^2);
%         rms = sum((GrF*u-Bdata).^2);
%         redu_perc = 100*(rms0-rms)/rms0;   
%         fprintf('rms misfit (dat., res.) = %e %e (%f%%) \n',rms0,rms,redu_perc);
%         fprintf('resnorm, resid. = %e %e \n',sqrt(resnorm),mean(residual));
        fprintf('exitflag is %d\n',exitflag);   % 1 means the function converged to a solution x   
        
        slip_model(:,12) = u(1:sum(tSm));
        slip_model(:,13) = u(sum(tSm)+1:end);
%         show_slip_model(slip_model);
              
        S1A_misfit = sum((S1A_Bdata - S1A_GrF * u).^2);
        test_misfit = sum((test_Bdata - test_GrF * u).^2);
        L_curve(ii,1) = S1A_misfit;
        L_curve(ii,2) = test_misfit;
        
    end
    save(['data_sensitivity/S1A_',test_data_type,'.mat'],'L_curve');
    figure; hold on
    plot(L_curve(:,1),L_curve(:,2),'r*-','linewidth',2);
    xlabel('\chi^2_{Sentinel-1 LOS}');
    ylabel(['\chi^2_{',test_data_type,' LOS}']);
    set(gca,'fontsize',20);
    set(gcf,'PaperPositionMode','auto');
end