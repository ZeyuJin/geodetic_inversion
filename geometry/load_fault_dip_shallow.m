function slip_model = load_fault_dip_shallow(fault_segment_file,varargin)
% The format of fault_segment_file should be like this:
% [xpt(1),ypt(1)]  [xpt(2),ypt(2)]  +  [xd(1),yd(1)]   [xd(2),yd(2)]
% surface trace   +  dip trace projected at the surface
% the default dipping depth is 20km, user could change it using varargin

% The function will discretize the dipping fault into multiple segments,
% even they have the same strike, but they are denoted with different
% segment ID, because they are not seamlessly adjacent with each other, 
% and it is convenient to add smoothing later on

% it is better to use 'dip_depth' and 'dip_change_id' at the same time,
% and they should have the same length (1 id : 1 depth)
% written by Zeyu Jin on Oct. 2019
   
   d2r = pi / 180;
   lon_eq = -117.5;
   lat_eq = 35.5;
   ref_lon = lon_eq;
   
   %% Default geometric values for the shallow dipping fault
   zstart = 0;
   dip_change_id = [];
   this_depth = 3e3;   % best fit for the northern/central part
   dip_depth = [];
   lp_top = 1e3;
   bias_lp = 1.0;
   bias_wp = 1.0;
   slip_model = [];
   fault_id = 0;
   
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
                       fault_id = varargin{CC*2};        % which segment start to count
                   case 'dip_change_id'
                       dip_change_id = varargin{CC*2};   % which dipping segment could be changed the dipping depth                       
                   case 'dip_depth'
                       dip_depth = varargin{CC*2};       
                   case 'len_top'
                       lp_top = varargin{CC*2};
                   case 'l_ratio'
                       bias_lp = varargin{CC*2};
                   case 'w_ratio'
                       bias_wp = varargin{CC*2};
                   case 'depth_start'
                       zstart = varargin{CC*2};
               end
           catch
               error('Unrecognized Keyword\n');
           end
       end
   end
  
   fault_data = load(fault_segment_file);
   nflt = size(fault_data,1);              % number of dipping faults   
%    [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);
   [xo,yo] = ll2xy(lon_eq,lat_eq,ref_lon);
   
   % read depth of different dipping segments
   if length(dip_depth) ~= length(dip_change_id) && length(dip_depth) ~= 1
       error('Not equal size between dipping segments and depth!');
   elseif max(dip_change_id) > nflt
       error('The NO.ID of dipping segment is out of range!');
   else
       tmp = dip_depth;
       dip_depth = this_depth .* ones(1,nflt);
       dip_depth(dip_change_id) = tmp;       
   end
   
   for ii = 1:nflt
       % load the UTM coordinates of this strike of dipping fault
       this_dip_segment = fault_data(ii,:);
       surface_lon_pt = [this_dip_segment(1);this_dip_segment(3)];
       surface_lat_pt = [this_dip_segment(2);this_dip_segment(4)];
       dip_lon_pt = [this_dip_segment(5);this_dip_segment(7)]; 
       dip_lat_pt = [this_dip_segment(6);this_dip_segment(8)];
%        [surface_xutm,surface_yutm] = utm2ll(surface_lon_pt,surface_lat_pt,0,1);
%        [dip_xutm,dip_yutm] = utm2ll(dip_lon_pt,dip_lat_pt,0,1);
       [surface_xutm,surface_yutm] = ll2xy(surface_lon_pt,surface_lat_pt,ref_lon);
       [dip_xutm,dip_yutm] = ll2xy(dip_lon_pt,dip_lat_pt,ref_lon);
       
       surface_xutm = surface_xutm - xo;  surface_yutm = surface_yutm - yo;
       dip_xutm = dip_xutm - xo;          dip_yutm = dip_yutm - yo;
       
       % the fault starts from the bottom to the top
       dx = surface_xutm(1) - surface_xutm(2); 
       dy = surface_yutm(1) - surface_yutm(2);
       length_this_segment = sqrt(dx^2 + dy^2);
       N_divide_segment = round(length_this_segment/lp_top);
       if N_divide_segment == 0, N_divide_segment = 1; end     % in case that this segment is too short

       L = length_this_segment / N_divide_segment;
       theta = atan2(dy,dx);
       strike = 90 - theta / d2r;
       if strike < 0, strike = strike + 360; end

       this_depth = dip_depth(ii);
%        disp(this_depth);
       
       for jj = 1:N_divide_segment      % compute the geometry of each segment with same strike direction
           fault_id = fault_id + 1;
           
           W_front = sqrt((surface_xutm(1)-dip_xutm(1))^2 + ...
                          (surface_yutm(1)-dip_yutm(1))^2);
           W_tail =  sqrt((surface_xutm(2)-dip_xutm(2))^2 + ... 
                          (surface_yutm(2)-dip_yutm(2))^2);
           W_proj =  (W_front-W_tail) / (2*N_divide_segment) * (2*jj-1) + W_tail;   % width of patch projected to the surface
           W = sqrt(W_proj^2 + this_depth^2);
           dip = asin(this_depth/W) / d2r;
           N_layer = floor(W/lp_top);
           if N_layer == 0, N_layer = 1; end
          
           xstart = (surface_xutm(1)-surface_xutm(2)) / N_divide_segment * (jj-1) + surface_xutm(2);
           ystart = (surface_yutm(1)-surface_yutm(2)) / N_divide_segment * (jj-1) + surface_yutm(2);
           model_segment = make_fault_segments(fault_id,xstart,ystart,zstart,strike,dip,L,W,N_layer,lp_top,bias_lp,bias_wp);
           slip_model = [slip_model;model_segment];
       end
   end
end