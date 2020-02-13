function [H,h1,indx_less_smooth_patch] = build_smooth_function(slip_model_vs,slip_model_ds,segment_file,intersect_file,ramp_choice,varargin)
% build the smooth function using first-order Tikhonov regularization
% First smooth within each fault segment
% Second smooth conjunction between fault segments
% Third smooth the intersection between the shallow dipping fault and vertical fault segment
% Option: apply less (or more) smoothing to some dipping features (standard as ratio 1)
   
   slip_model = [slip_model_vs;slip_model_ds];
   slip_model(:,2)=[1:size(slip_model,1)]';   % recomputed finally to combine all the fault segments
   
   [H_plane,~] = smoo1_each_plane(slip_model);
   if ~isempty(segment_file)
       H_segment = smoo1_segments(slip_model,segment_file);
   else
       H_segment = [];
   end
   if ~isempty(intersect_file)
      H_intersect = smoo1_dip_vertical_intersection(slip_model_ds,slip_model_vs,intersect_file);
   else
      H_intersect = [];
   end
   H = [H_plane;H_segment;H_intersect];
   h1 = size(H,1);
   Np = size(H,2) / 2;   % number of patches
   
   dip_segments_id = [];
   % specify which dipping segments to be less smoothing for the dip component
   if ~isempty(varargin)
       for CC = 1:floor(length(varargin)/2)
           try
               switch lower(varargin{CC*2-1})
                   case 'dip_id'
                       dip_segments_id = varargin{CC*2};
               end
           catch
               error('Unrecognized Keyword');
           end
       end
   end
   
   % apply 1 to the corresponding dipping components same as strike-slip
   % ii -- loop through each row of smooth matrix H
   % jj -- loop through each dipping segments
   % kk -- loop through each patch
   indx_less_smooth_patch = [];
   if ~isempty(dip_segments_id)
      len_dip_seg = length(dip_segments_id);
      all_fault_id = slip_model(:,1);
      all_patch_id = slip_model(:,2);
%       all_layer_id = slip_model(:,3);
      
      indx_dip_patch = [];
      for jj = 1:len_dip_seg
      % find corresponding patches of shallow dip fault
          this_dip_segment_id = dip_segments_id(jj);
          indx_dip_segment = all_fault_id == this_dip_segment_id;
          tmp_patch = all_patch_id(indx_dip_segment);
%           tmp_layer = all_layer_id(indx_dip_segment);
%           indx_shallow_two = tmp_layer <= 2;      % top two layers          
%           tmp = tmp_patch(indx_shallow_two);
          indx_dip_patch = [indx_dip_patch;tmp_patch];
      end
%       disp(indx_dip_patch);
%       indx_vertical_patch = [1:3,21:23,36:38,76:80,87:88]';
        indx_vertical_patch = [1:3,21:23,36:38,76:80]';
%       indx_vertical_patch = [1:3,21:23,76:79]';
        indx_less_smooth_patch = [indx_vertical_patch;indx_dip_patch];
        tmp = [85:111];
        indx_more_smooth_patch = [tmp,tmp+Np]';    % add more smoothing of segment 3
      
      for ii = 1:h1
          indx_nonzero = find(H(ii,:));
          if length(indx_nonzero) ~= 2
             error('There are more than 2 elements in each row of smoothing matrix');
          end
          
          % reduce the partition effect on segment 3
          if ismember(indx_nonzero(1),indx_more_smooth_patch) || ismember(indx_nonzero(2),indx_more_smooth_patch)
              H(ii,indx_nonzero) = H(ii,indx_nonzero) .* 2;
          end          
          
          if max(indx_nonzero) < Np, continue; end    % skip the constraint on the strike component          
          left_nonzero = indx_nonzero(1) - Np;        % find only dip components
          right_nonzero = indx_nonzero(2) - Np;  
          if ismember(left_nonzero,indx_less_smooth_patch) && ismember(right_nonzero,indx_less_smooth_patch)
              ratio = max(abs(H(ii,indx_nonzero)));
              H(ii,indx_nonzero) = H(ii,indx_nonzero) ./ 5; % ratio;
%               disp(indx_nonzero-Np);
          end
          
%           % check the smoothing matrix of modelB
%           if indx_nonzero(1)-Np == 449 || indx_nonzero(2)-Np == 449 || indx_nonzero(1)-Np == 450 || indx_nonzero(2)-Np == 450
%               disp(indx_nonzero);
%           end
      end
   end  
      
   ramp_choice = lower(ramp_choice);   
   if strcmp(ramp_choice,'bi_ramp')
       rmp = zeros(h1,4);       % bilinear ramp
   elseif strcmp(ramp_choice,'qu_ramp_7')
       rmp = zeros(h1,7);
   elseif strcmp(ramp_choice,'qu_ramp_5')
       rmp = zeros(h1,5);
   else
       rmp = [];
   end
   H = [H,rmp];
   
%    % just for test smoothing index
%    for kk = 1:h1
%        tmp = find(H(kk,:));
%        if tmp(1) == 81 || tmp(2) == 81
%            disp(tmp);
%        end
%    end
end