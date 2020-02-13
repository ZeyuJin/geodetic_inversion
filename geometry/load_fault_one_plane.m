function slip_model = load_fault_one_plane(fault_segment_file,varargin)
% this function is used to build up multiple whole fault plane
% the fault plane can be vertical or dipping
% could not be used for shallow irregular dipping feature

    d2r=pi/180;
    lon_eq = -117.5;
    lat_eq = 35.5;

    %% Default geometric values for the vertical fault
    slip_model = [];
    fault_id = 0;
    % divide fault patch versus depth
    W=25e3;       % make the faults deeper  
    N_layer=8;    
    lp_top=1e3;   % should be the similar spatial size with sampled data
    bias_lp=1.3;
    bias_wp=1.3;
    zstart=0;
    this_dip = 90;
    dip = [];
    dip_change_id = [];
    
    %% read varargin values and assembly
    if ~isempty(varargin)
        for CC = 1:floor(length(varargin)/2)
            try
                switch lower(varargin{CC*2-1})
                    case 'lonc'
                        lon_eq = varargin{CC*2};
                    case 'latc'
                        lat_eq = varargin{CC*2};
                    case 'fault_id'
                        fault_id = varargin{CC*2};    % which segment start to count
                    case 'dip_change_id'
                        dip_change_id = varargin{CC*2}; 
                    case 'width'
                        W = varargin{CC*2};
                    case 'layers'
                        N_layer = varargin{CC*2};
                    case 'len_top'
                        lp_top = varargin{CC*2};
                    case 'l_ratio'
                        bias_lp = varargin{CC*2};
                    case 'w_ratio'
                        bias_wp = varargin{CC*2};
                    case 'depth_start'
                        zstart = varargin{CC*2};
                    case 'dip'
                        dip = varargin{CC*2};
                end
            catch
                error('Unrecognized Keyword\n');
            end
        end
    end
    
    %% read fault data
    [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);
    fault_data = load(fault_segment_file);
    lon_pt = [fault_data(:,1);fault_data(:,3)];
    lat_pt = [fault_data(:,2);fault_data(:,4)];
    nflt = size(fault_data,1);

    [xutm_pt,yutm_pt]=utm2ll(lon_pt,lat_pt,0,1);
    xpt=xutm_pt-xo;
    ypt=yutm_pt-yo;   
    
    % change dip angle for each segment
    if length(dip) ~= length(dip_change_id) && length(dip) ~= 1
        error('Not equal size between dipping segments and angle!');
    elseif max(dip_change_id) > nflt
        error('The NO.ID of dipping segment is out of range!');
    else
        tmp = dip;
        dip = this_dip .* ones(1,nflt);
        dip(dip_change_id) = tmp;
    end
    
    %% the fault geometry is defined as Peter Shearer's textbook   
    for kk=1:nflt
        xstart = xpt(kk+nflt);
        ystart = ypt(kk+nflt);   
        xend = xpt(kk);
        yend = ypt(kk);   
        dx=xend-xstart;      % negative constraint (the fault starts from the bottom to the top)
        dy=yend-ystart;   
   
        L=sqrt(dx^2+dy^2);
        theta=atan2(dy,dx);
        strike1=90-theta/d2r;
        dip_angle = dip(kk);
%         disp(dip_angle);
   
        if (strike1<0)
            strike1=strike1+360;
        end

        fault_id = fault_id + 1;
        model_segment = make_fault_segments(fault_id,xstart,ystart,zstart,strike1,dip_angle,L,W,N_layer,lp_top,bias_lp,bias_wp);
        slip_model=[slip_model;model_segment];
    end
end