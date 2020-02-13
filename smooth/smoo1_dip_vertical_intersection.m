function smooth = smoo1_dip_vertical_intersection(slip_model_ds,slip_model_vs,connection_file,varargin)
% apply 1st derivative smoothing at the intersection between dip and vertical faults
% user should provide the segment ID of both dipping and vertical faults in one file
% to find the connection between shallow dip part and the vertical segment
% e.g.      1  5
%           2  5
%           3  5
%           4  9
% the first column means the segment_ID of dipping faults
% the second column means the segment_ID of vertical faults
% solve the bug for the uniform patches on Nov. 2019

    connect_id = load(connection_file);
    smooth = [];
    slip_model = [slip_model_vs;slip_model_ds];
    Np = size(slip_model,1);      % number of total patches
    Nv = size(slip_model_vs,1);   % number of vertical patches
    Fdip = 5;
    if ~isempty(varargin)
       for CC = 1:floor(length(varargin)/2)
           try
               switch lower(varargin{CC*2-1})
                   case 'fdip'
                       Fdip = varargin{CC*2};
               end
           catch
               error('Unrecognized Keyword');
           end
       end
    end
    
    all_strike = slip_model_ds(:,9);
    strike_uniq = unique(all_strike,'stable');  % reduce the redundancy of strike information
    % also same order with index of dipping segments ID, e.g. [1 2 3 4 ...]
        
    % the following loop index ii, jj mean:
    % ii -- which pair of connection (dip + vertical)
    % jj -- which small dipping part of this segment
    % mm -- which patch in the bottom layer of this dipping part
    % kk -- which patch in one layer of vertical segment to be intersected with each dipping patch
    
    for ii = 1:size(connect_id,1)
        dip_con_id = connect_id(ii,1);
        stk_con_id = connect_id(ii,2);
        
        strike_this_dip_segment = strike_uniq(dip_con_id);
        indx_dip_segment = find(all_strike == strike_this_dip_segment);  
        % use the strike information of dipping faults to find all small segments
        dip_this_segment = slip_model_ds(indx_dip_segment,:);
        dip_ID_this_segment = dip_this_segment(:,1);       % the real ID of each small dipping segments
        dip_ID_uniq_this_segment = unique(dip_ID_this_segment,'stable');   % reduce the redundancy of dip segment ID
        
        indx_stk_segment = find(slip_model_vs(:,1) == stk_con_id);
        vertical_this_segment = slip_model_vs(indx_stk_segment,:);
        
        for jj = 1:length(dip_ID_uniq_this_segment)
            this_dip_ID = dip_ID_uniq_this_segment(jj);
            indx_all_patch_this_fault = find(dip_ID_this_segment == this_dip_ID);
            
            % find the bottom patch of this small dipping fault
            dip_this_fault = dip_this_segment(indx_all_patch_this_fault,:);
            bottom_layer = max(dip_this_fault(:,3));
%             disp(bottom_layer);
            W_all_patch = zeros(1,bottom_layer);
            dip_angle_all_patch = zeros(1,bottom_layer);
            for nn = 1:bottom_layer
                tmp = find(dip_this_fault(:,3)==nn,1);     % find the first patch of bottom layer
                W_all_patch(nn) = dip_this_fault(tmp,8);
                dip_angle_all_patch(nn) = dip_this_fault(tmp,10);
            end
            intersect_depth = (-1) * sum(W_all_patch .* sind(dip_angle_all_patch));    % the dipping depth of this fault
            
            % loop for each bottom patch of dipping fault
            N_patch_bottom_layer = length(find(dip_this_fault(:,3) == bottom_layer));
            indx_patch_bottom_layer = find(dip_this_fault(:,3) == bottom_layer);
            for mm = 1:N_patch_bottom_layer
                indx_this_bottom_patch = indx_patch_bottom_layer(mm);
                bottom_patch = dip_this_fault(indx_this_bottom_patch,:);   % the bottom patch (the indx end may be improved)
                xp_bottom_patch = bottom_patch(4);
                yp_bottom_patch = bottom_patch(5);
                zp_bottom_patch = bottom_patch(6);
                lp_bottom_patch = bottom_patch(7);
                wp_bottom_patch = bottom_patch(8);
                stk_bottom_patch = strike_this_dip_segment;
                dip_bottom_patch = bottom_patch(10);
            
                % compute the center coordinates pf the bottom patch
                [xc_dip,yc_dip,zc_dip] = compute_patch_coords(xp_bottom_patch,yp_bottom_patch,zp_bottom_patch, ...
                    lp_bottom_patch,wp_bottom_patch,stk_bottom_patch,dip_bottom_patch);
                       
                % find the intersection layer of this fault segment
                zp_vertical = unique(vertical_this_segment(:,6),'stable');
                wp_vertical = unique(vertical_this_segment(:,8),'stable');    % has only one value for uniform patches
                dip_vertical = unique(vertical_this_segment(:,10),'stable');  % const for one plane fault
                if length(wp_vertical) == 1, wp_vertical = wp_vertical .* ones(length(zp_vertical),1); end
                num_of_layers = max(vertical_this_segment(:,3));
%                 disp(zp_vertical'/1000);
%                 disp(wp_vertical'/1000);
                for count = 1:num_of_layers
                    z_top_this_layer =  zp_vertical(count);
                    z_bottom_this_layer = zp_vertical(count) - wp_vertical(count) .* sind(dip_vertical);         
                    if (z_top_this_layer > intersect_depth) && (intersect_depth >= z_bottom_this_layer)
%                         disp([z_top_this_layer/1000  z_bottom_this_layer/1000]);
%                         disp(['Intersect with the ',num2str(count),' layer in No.',num2str(stk_con_id), ' vertical segment.']);
                        indx_this_layer = find(vertical_this_segment(:,3) == count);
                        vertical_this_layer = vertical_this_segment(indx_this_layer,:);
                        break
                    end
                end
            
                % compute the distance between bottom patch in dipping fault
                % and adjacent patch in vertical fault
                num_of_patches_this_layer = length(indx_this_layer);
                dis_dip_stk = zeros(1,num_of_patches_this_layer);  % save the distance between two patches
                for kk = 1:num_of_patches_this_layer
                    xp_this_patch = vertical_this_layer(kk,4);
                    yp_this_patch = vertical_this_layer(kk,5);
                    zp_this_patch = vertical_this_layer(kk,6);
                    lp_this_patch = vertical_this_layer(kk,7);
                    wp_this_patch = vertical_this_layer(kk,8);
                    stk_this_patch = vertical_this_layer(kk,9);
                    dip_this_patch = vertical_this_layer(kk,10);
                
                    [xc_stk,yc_stk,zc_stk] = compute_patch_coords(xp_this_patch,yp_this_patch,zp_this_patch, ...
                        lp_this_patch,wp_this_patch,stk_this_patch,dip_this_patch);
                    dis_dip_stk(kk) = p2p_dis([xc_dip,yc_dip,zc_dip],[xc_stk,yc_stk,zc_stk]);
                end
                % find the closest distance and compute the index of two
                % patches in both slip_model_vs || slip_model_ds
                col = zeros(2,Np*2);
               
                [~,n1] = sort(dis_dip_stk);
                indx_this_patch = n1(1);         %  the first one is the closest one
                indx_above_this_layer = length(find(vertical_this_segment(:,3) < count));
                indx_pre_this_segment = length(find(slip_model_vs(:,1) < stk_con_id));
                indx_patch_vertical = indx_this_patch + indx_above_this_layer + indx_pre_this_segment;   % index in slip_model_vs
            
                indx_dip_patch = mm;
                indx_above_dip_layer = length(find(dip_this_fault(:,3) < bottom_layer));
                indx_pre_this_fault = length(find(slip_model_ds(:,1) < this_dip_ID));
                indx_patch_dip = indx_dip_patch + indx_above_dip_layer + indx_pre_this_fault + Nv;       % index in slip_model_ds
%                 disp([indx_dip_patch,indx_above_dip_layer,indx_pre_this_fault]);
                
                col(1,indx_patch_vertical) = -1;
                col(1,indx_patch_dip) = 1;
                col(2,indx_patch_vertical+Np) = -Fdip;
                col(2,indx_patch_dip+Np) = Fdip * sind(dip_bottom_patch);    % Is it right?
%                 disp([indx_patch_vertical,indx_patch_dip,indx_patch_vertical+Np,indx_patch_dip+Np]);
                smooth = [smooth;col];
            end
        end
    end
end