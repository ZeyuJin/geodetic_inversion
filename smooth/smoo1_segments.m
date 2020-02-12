function smooth = smoo1_segments(slip_model_in,segment_file,varargin)
% apply the 1st derivative smoothing
% user should provide the segment file to clarify that which segments are
% connected with each other
% the segment file could be in the format that:
% (1 right 2 left) | (1 right 2 right) | (1 left 2 left)

   fault_id = slip_model_in(:,1);
   indx_layer = slip_model_in(:,3);
   wp = slip_model_in(:,8);
   dip = slip_model_in(:,10);
   
   Np = length(indx_layer);                 % number of patches
   Nf = max(fault_id);                      % number of fault segments
   Nl = max(indx_layer);                    % number of layers 
   WP = zeros(1,Nf);                        % number of layers of each segment
   for ii = 1:Nf
       indx_this_segment = find(fault_id == ii);
       layer_this_segment = indx_layer(indx_this_segment);
       WP(ii) = max(layer_this_segment);
   end
   
   nL = zeros(Nf,Nl);                       % number of patches in each segment/each layer
   dW = zeros(Nf,Nl);                       % width of patches in each segment/each layer
   
   for ii = 1:Nf
       indx_this_segment = find(fault_id == ii);
       layer_this_segment = indx_layer(indx_this_segment);
       wp_this_segment = wp(indx_this_segment);
       for jj = 1:WP(ii)
           tmp = find(layer_this_segment == jj);
           nL(ii,jj) = length(tmp);
           Wpatch = wp_this_segment(tmp);
           dW(ii,jj) = Wpatch(1);
       end
   end
   
   tSm = sum(nL,2)';
   tSm = [0,tSm];    % number of patches of each segment [1x(Nf+1)]
   Fdip = 5;
   smooth = [];
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
   
   % open the segment_file to find the two adjacent segments
   fid = fopen(segment_file);
   C = textscan(fid,'%d %s %d %s');
   seg_id1 = C{:,1};   seg_id2 = C{:,3};
   side1 = C{:,2};     side2 = C{:,4};
   fclose(fid);
   
   for count = 1:length(seg_id1)
       id1 = seg_id1(count);
       id2 = seg_id2(count);
       s1 = side1(count);
       s2 = side2(count);
       
       dip_id1 = dip(find(fault_id == id1,1));
       dip_id2 = dip(find(fault_id == id2,1));
       
       for jj_first = 1:WP(id1)
           if strcmp(s1,'right')
               indx_minus = sum(tSm(1:id1)) + sum(nL(id1,1:jj_first));
           else
               indx_minus = sum(tSm(1:id1)) + sum(nL(id1,1:jj_first-1)) + 1;
           end
           wp_first_top = sum(dW(id1,1:jj_first-1));   % * sind(dip_id1);  % use the depth not width to find the contact
           wp_first_down = sum(dW(id1,1:jj_first));    % * sind(dip_id1);
           
           for jj_second = 1:WP(id2)
               if strcmp(s2,'right')
                   indx_plus = sum(tSm(1:id2)) + sum(nL(id2,1:jj_second));
               else
                   indx_plus = sum(tSm(1:id2)) + sum(nL(id2,1:jj_second-1)) + 1;
               end
               wp_second_top = sum(dW(id2,1:jj_second-1));   % * sind(dip_id2);
               wp_second_down = sum(dW(id2,1:jj_second));    % * sind(dip_id2);
               
               % need to be changed during every iteration to find the contact, because there maybe multiple contacts
               col = zeros(2,2*Np);    
               if (abs(wp_first_top - wp_second_top) < 1e-3) || (abs(wp_first_down - wp_second_down) < 1e-3) ...
                    || (wp_first_top < wp_second_top && wp_second_top < wp_first_down) ...
                    || (wp_first_top < wp_second_down && wp_second_down < wp_first_down)   % has contact
                    indx_minus_strike = indx_minus;
                    indx_minus_dip = indx_minus_strike + Np;
                    indx_plus_strike = indx_plus;
                    indx_plus_dip = indx_plus_strike + Np;
                    col(1,indx_minus_strike) = -1;
                    col(2,indx_minus_dip) = -Fdip * sind(dip_id1);
                    col(1,indx_plus_strike) = 1;
                    col(2,indx_plus_dip) = Fdip * sind(dip_id2);                    
                    smooth = [smooth;col];
%                     disp([indx_minus_strike,indx_plus_strike,indx_minus_dip,indx_plus_dip]);
               end
           end
       end
   end
end