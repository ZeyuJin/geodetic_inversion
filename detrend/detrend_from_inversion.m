function detrend_from_inversion(slip_model_vs,slip_model_ds,detrend_file,varargin)

    %% open the detrend list of InSAR data
    fid = fopen(detrend_file);       % detrend_file = 'detrend.list';
    C = textscan(fid,'%s %d\n');
    track = C{:,1};
    ntrack = length(track);
    fclose(fid);

    %% the postivity constraint (adapted from Yuri's code)
    slip_model = [slip_model_vs;slip_model_ds];
    slip_model(:,2)=[1:size(slip_model,1)]';    % recomputed finally to combine all the fault segments
    nflt = max(slip_model(:,1));
    tSm = zeros(1,nflt+1);
    fault_id = slip_model(:,1);
    for ii=1:nflt
        tSm(ii+1) = length(find(fault_id == ii));
    end
    NT = 2; NS = nflt;  % the number of segments  
%     fault_patch = slip_model(:,1:3);

    %% some parameters in the inversion
    beta = 5;                               % relative weighting between GPS and InSAR as a default
    lambda = 8e-2;                          % weight of smoothness factor
    ramp_array = 4 * ones(1,ntrack);        % array of ramp index number, should be the same size with detrend.list
    segment_file = [];      intersect_file = [];      % segment_file = 'seg_connect_coarse';    intersect_file = 'dip_vertical_connect';
    data_type = 'insar';
    model_type = 'okada';
    
    %% read varargin values and assembly
    if ~isempty(varargin)
        for CC = 1:floor(length(varargin)/2)
            try
                switch lower(varargin{CC*2-1})
                    case 'gps2insar'
                        beta = varargin{CC*2};
                    case 'smoothness'
                        lambda = varargin{CC*2};
                    case 'ramp_order'
                        ramp_array = varargin{CC*2};
                        if size(ramp_array,1) ~= 1 || size(ramp_array,2) ~= ntrack
                            error('Ramp Order not correct!');
                        end
                    case 'segment_smooth_file'
                        segment_file = varargin{CC*2};
                    case 'intersect_smooth_file'
                        intersect_file = varargin{CC*2};
                    case 'data_type'
                        data_type = varargin{CC*2};
                    case 'model_type'
                        model_type = varargin{CC*2};    % homogenous or layered    
                end
            catch
                error('Unrecognized Keyword');
            end
        end
    end


    %% each interferogram/AZO has to be de-trended individually
    for kk = 1:length(ramp_array)
        this_track = track{kk};
        [ramp_type,add_col] = detrend_ramp_type(ramp_array(kk));
        [Gl_raw,Gl,bdl_raw,bdl] = build_green_function(slip_model,[this_track,'/los_samp_uniform.mat'],data_type,ramp_type,model_type);
        [Gp_raw,Gp,bp_raw,bp] = build_green_function(slip_model,'GPS/Floyd_GPS/Floyd_GPS_all.mat','camp_gps',ramp_type,model_type,beta);
        [Gs_raw,Gs,bs_raw,bs] = build_green_function(slip_model,'GPS/survey_GPS/survey_gps_2d.mat','camp_gps',ramp_type,model_type,3);
        [H,h1] = build_smooth_function(slip_model_vs,slip_model_ds,segment_file,intersect_file,ramp_type);
    
        G_raw = [Gl;Gp;Gs];                    Greens = [G_raw;H*lambda/h1];
        bd_raw = [bdl;bp;bs];                  bdata_sm = [bd_raw;zeros(h1,1)];
        [lb,ub] = bounds_coarse_A(NS,NT,tSm,add_col,slip_model);   % do not add any boundary condition

        % linear inversion
        options = optimset('LargeScale','on','DiffMaxChange',1e-1,'DiffMinChange',1e-12, ...
            'TolCon',1e-12,'TolFun',1e-12,'TolPCG',1e-12,'TolX',1e-12,'MaxIter',1e9,'MaxPCGIter',1e9);
        [u,resnorm,residual,exitflag] = lsqlin(Greens,double(bdata_sm),[],[],[],[],lb,ub,[],options);
    
        % save in the slip_model
        slip_model(:,12) = u(1:sum(tSm));
        slip_model(:,13) = u(sum(tSm)+1:2*sum(tSm));
        ramp_coef = u(end-add_col+1:end);

        %% plot and save the finite fault inversion
%         slip_model(:,2)=[1:size(slip_model,1)]';
        show_slip_model(slip_model);
        write_slip_model_okada(slip_model,'fault_initial_coarse.slip');

        %% plot the resampled data fitting
%         insar_model = Gl_raw * u;
        cgps_model = Gp_raw * u;
        survey_gps_model = Gs_raw * u;

%         plot_insar_model_resampled([this_track,'/los_samp_uniform.mat'],insar_model);
        plot_gps_model('GPS/Floyd_GPS/Floyd_GPS_all.mat',cgps_model,'data_type','survey');
        plot_gps_model('GPS/survey_GPS/survey_gps_2d.mat',survey_gps_model,'data_type','survey','site_name','GPS/survey_GPS/campaign_sites');

        %% plot the original/de-trended interferograms and ramp
        plot_insar_detrend(this_track,'los_ll_uniform.grd','dem_uniform.grd',ramp_coef);
    end

end