function [slip_model,fault_id] = load_fault_dip(fault_segment_file,fault_id)
    addpath('/Users/zej011/work/Kang_tutorial/codes_utilities/matlab/igppsar');
    
    d2r = pi/180;
    lon_eq = -117.5;
    lat_eq = 35.5;
    [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);
    
    fault_data = load(fault_segment_file);
    lon_pt = [fault_data(:,1);fault_data(:,3)];
    lat_pt = [fault_data(:,2);fault_data(:,4)];
    nflt = size(fault_data,1);
    
    [xutm_pt,yutm_pt]=utm2ll(lon_pt,lat_pt,0,1);
    xpt=xutm_pt-xo;
    ypt=yutm_pt-yo;
    
    % the coordinates of connected deep fault
    lon1_start = -117.619167;   lat1_start = 35.797283;
    lon1_end = -117.444167;     lat1_end = 35.634303;
    lon1_mid1 = -117.544501;    lat1_mid1 = 35.727822;
    lon1_mid2 = -117.468756;    lat1_mid2 = 35.657349;
    
    [x1_start,y1_start] = utm2ll(lon1_start,lat1_start,0,1);
    x1_start = x1_start - xo;  y1_start = y1_start - yo;
    
    [x1_end,y1_end] = utm2ll(lon1_end,lat1_end,0,1);
    x1_end = x1_end - xo;      y1_end = y1_end - yo;
    
    [x1_mid1,y1_mid1] = utm2ll(lon1_mid1,lat1_mid1,0,1);
    x1_mid1 = x1_mid1 - xo;    y1_mid1 = y1_mid1 - yo;
    
    [x1_mid2,y1_mid2] = utm2ll(lon1_mid2,lat1_mid2,0,1);
    x1_mid2 = x1_mid2 - xo;    y1_mid2 = y1_mid2 - yo;
    
%    %% first dipping segment
    slip_model = [];
    
    zstart = 0;
    dip_depth = 2e3;
    lp_top = 1e3;
    bias_lp = 1.0;
    bias_wp = 1.0;
    
    fault_id = fault_id + 1;
    dx = xpt(1)-xpt(1+nflt);     % positive constraint
    dy = ypt(1)-ypt(1+nflt);
    this_segment = sqrt(dx^2 + dy^2);
    N_segment = round(this_segment / lp_top);
    L = this_segment / N_segment;
    theta = atan2(dy,dx);
    strike = 90 - theta / d2r;
    if (strike < 0)
        strike = strike + 360;
    end
    
    for ii = 1:N_segment        
        W_max = sqrt((xpt(2)-x1_mid1)^2 + (ypt(2) - y1_mid1)^2);
        W_proj = W_max / 2 / N_segment * (2*(N_segment+1-ii)-1);     % horizontal surface length
        W = sqrt(W_proj^2 + dip_depth^2); 
        dip = asin(dip_depth / W) / d2r;
        N_layer = floor(W/lp_top);
        if N_layer == 0
            N_layer = 1;
        end
        
        xstart = (xpt(1) - xpt(1+nflt)) / N_segment * (ii-1) + xpt(1+nflt);
        ystart = (ypt(1) - ypt(1+nflt)) / N_segment * (ii-1) + ypt(1+nflt);
        model_segment = make_fault_segments(fault_id,xstart,ystart,zstart,strike,dip,L,W,N_layer,lp_top,bias_lp,bias_wp);
        slip_model=[slip_model;model_segment];
    end
    
%    %% second dipping segment
    fault_id = fault_id + 1;
    dx = xpt(2)-xpt(2+nflt);     % positive constraint
    dy = ypt(2)-ypt(2+nflt);
    this_segment = sqrt(dx^2 + dy^2);
    N_segment = round(this_segment / lp_top);
    L = this_segment / N_segment;
    theta = atan2(dy,dx);
    strike = 90 - theta / d2r;
    if (strike < 0)
        strike = strike + 360;
    end
    
    for ii = 1:N_segment       
        W_max = sqrt((xpt(2)-x1_mid1)^2 + (ypt(2) - y1_mid1)^2);
        W_min = sqrt((xpt(3)-x1_mid2)^2 + (ypt(3) - y1_mid2)^2);
        W_proj = (W_max - W_min) / 2 / N_segment * (2*ii-1) + W_min;    % horizontal surface length
        W = sqrt(W_proj^2 + dip_depth^2);
        dip = asin(dip_depth / W) / d2r;
        N_layer = floor(W/lp_top);
        if N_layer == 0
            N_layer = 1;
        end
    
        xstart = (xpt(2) - xpt(2+nflt)) / N_segment * (ii-1) + xpt(2+nflt);
        ystart = (ypt(2) - ypt(2+nflt)) / N_segment * (ii-1) + ypt(2+nflt);
        model_segment = make_fault_segments(fault_id,xstart,ystart,zstart,strike,dip,L,W,N_layer,lp_top,bias_lp,bias_wp);
        slip_model = [slip_model;model_segment];
    end

%    %% third dipping segment
    fault_id = fault_id + 1;
    dx = xpt(3) - xpt(3+nflt);
    dy = ypt(3) - ypt(3+nflt);
    this_segment = sqrt(dx^2 + dy^2);
    N_segment = round(this_segment / lp_top);
    L = this_segment / N_segment;
    theta = atan2(dy,dx);
    strike = 90 - theta / d2r;
    if (strike < 0)
        strike = strike + 360;
    end
    
    for ii = 1:N_segment
        W_max = sqrt((xpt(3)-x1_mid2)^2 + (ypt(3)-y1_mid2)^2);
        W_proj = W_max / 2 / N_segment * (2*ii-1);     % horizontal surface length
        W = sqrt(W_proj^2 + dip_depth^2);
        dip = asin(dip_depth / W) / d2r;
        N_layer = floor(W/lp_top);
        if N_layer == 0 
            N_layer = 1;
        end
        
        xstart = (xpt(3) - xpt(3+nflt)) / N_segment * (ii-1) + xpt(3+nflt);
        ystart = (ypt(3) - ypt(3+nflt)) / N_segment * (ii-1) + ypt(3+nflt);
        model_segment = make_fault_segments(fault_id,xstart,ystart,zstart,strike,dip,L,W,N_layer,lp_top,bias_lp,bias_wp);
        slip_model = [slip_model;model_segment];
    end
    
%    %% forth dipping segment 
     fault_id = fault_id + 1;
     lon2_start = -117.418725;  lat2_start = 35.611650;
     lon2_end = -117.223206;    lat2_end = 35.492238;
     
     [x2_start,y2_start] = utm2ll(lon2_start,lat2_start,0,1);
     x2_start = x2_start - xo;  y2_start = y2_start - yo;
     
     [x2_end,y2_end] = utm2ll(lon2_end,lat2_end,0,1);
     x2_end = x2_end - xo;      y2_end = y2_end - yo;
     
     dx = xpt(4) - xpt(4+nflt);
     dy = ypt(4) - ypt(4+nflt);
     this_segment = sqrt(dx^2 + dy^2);
     N_segment = round(this_segment / lp_top);
     L = this_segment / N_segment;
     theta = atan2(dy,dx);
     strike = 90 - theta / d2r;
     if (strike < 0)
         strike = strike + 360;
     end
     
     for ii = 1:N_segment
         W_max = sqrt((xpt(4+nflt)-x2_end)^2 + (ypt(4+nflt)-y2_end)^2);
         W_proj = W_max / 2 / N_segment * (2*(N_segment+1-ii)-1);   % horizontal surface length
         W = sqrt(W_proj^2 + dip_depth^2);
         dip = asin(dip_depth / W) / d2r;
         N_layer = floor(W/lp_top);
         if N_layer == 0 
             N_layer = 1;
         end
         
         xstart = (xpt(4) - xpt(4+nflt)) / N_segment * (ii-1) + xpt(4+nflt);
         ystart = (ypt(4) - ypt(4+nflt)) / N_segment * (ii-1) + ypt(4+nflt);
         model_segment = make_fault_segments(fault_id,xstart,ystart,zstart,strike,dip,L,W,N_layer,lp_top,bias_lp,bias_wp);
         slip_model = [slip_model;model_segment];
     end
end