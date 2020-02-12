function [smooth,nL] = smoo1_each_plane(slip_model_in,varargin)
% apply the 1st derivative smoothing
% assume each fault plane has the same number of layers
% solving the comparison between two floating numbers on Oct. 23rd, 2019
% also solving the index bug of nL (Nf -> ii) loop
% solving the smoothing index between strike and dip components
% col = zeros(2,Np*2); not to be only one row 

   fault_id = slip_model_in(:,1);
   indx_layer = slip_model_in(:,3);
   lp = slip_model_in(:,7);
   
   Np = length(indx_layer);                 % number of patches
   Nf = max(fault_id);                      % number of fault segments
   Nl = max(indx_layer);                    % number of layers 
   Wp = zeros(1,Nf);                        % number of layers of each segment
   for ii = 1:Nf
       indx_this_segment = fault_id == ii;
       layer_this_segment = indx_layer(indx_this_segment);
       Wp(ii) = max(layer_this_segment);
   end
   
   nL = zeros(Nf,Nl);                       % number of patches in each segment/each layer
   dL = zeros(Nf,Nl);                       % length of patches in each segment/each layer
   
   for ii = 1:Nf
       indx_this_segment = find(fault_id == ii);
       layer_this_segment = indx_layer(indx_this_segment);
       lp_this_segment = lp(indx_this_segment);
       for jj = 1:Wp(ii)
           tmp = find(layer_this_segment == jj);
           nL(ii,jj) = length(tmp);
           Lpatch = lp_this_segment(tmp);
           dL(ii,jj) = Lpatch(1);
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
   
   
   % the following loop index ii, jj, kk mean:
   % ii -- which fault segment
   % jj -- which layer of this segment
   % kk -- which patch in this layer of this segment
   
   % smooth between patches in the same layer and same fault segment
   for ii = 1:Nf  
       for jj = 1:Wp(ii)
           for kk = 1:nL(ii,jj)-1         % impose 1st derivative in length
               col = zeros(2,Np*2);       % strike + dip (2 rows)
               % patches in previous segments + patches at top layers of
               % this segment + number of this patch
               indx_minus_strike = sum(tSm(1:ii)) + sum(nL(ii,1:jj-1)) + kk;   
               indx_minus_dip = indx_minus_strike + Np;
               indx_plus_strike = sum(tSm(1:ii)) + sum(nL(ii,1:jj-1)) + (kk+1);
               indx_plus_dip = indx_plus_strike + Np;
               col(1,indx_minus_strike) = -1;
               col(2,indx_minus_dip) = -Fdip;
               col(1,indx_plus_strike) = 1;
               col(2,indx_plus_dip) = Fdip;
               smooth = [smooth;col];
%                disp([indx_minus_strike,indx_plus_strike,indx_minus_dip,indx_plus_dip]);
           end
       end
   end
   
   % smooth between patches in the adjacent layer and same fault segment
   for ii = 1:Nf
       for jj = 1:Wp(ii)-1     % from top layer to the (bottom-1) layer
           for kk = 1:nL(ii,jj)      
               
               Lp_top_left = (kk-1)*dL(ii,jj);
               Lp_top_right = kk*dL(ii,jj);
               for kk_down = 1:nL(ii,jj+1)
                   Lp_bottom_left = (kk_down-1)*dL(ii,jj+1);
                   Lp_bottom_right = kk_down*dL(ii,jj+1);
                   
                   % be careful with the "=" sign! 
                   % only left = left; right = right
                   col = zeros(2,Np*2);
                   if (abs(Lp_top_left - Lp_bottom_left) < 1e-3) || (abs(Lp_top_right - Lp_bottom_right) < 1e-3) ...   
                           || (Lp_top_left < Lp_bottom_left && Lp_bottom_left < Lp_top_right) ...         % can not only be adjacent at the boundary
                           || (Lp_top_left < Lp_bottom_right && Lp_bottom_right < Lp_top_right)           % has contact
                       indx_minus_strike = sum(tSm(1:ii)) + sum(nL(ii,1:jj-1)) + kk;
                       indx_minus_dip = indx_minus_strike + Np;
                       indx_plus_strike = sum(tSm(1:ii)) + sum(nL(ii,1:jj)) + kk_down;
                       indx_plus_dip = indx_plus_strike + Np;
                       col(1,indx_minus_strike) = -1;
                       col(2,indx_minus_dip) = -Fdip;
                       col(1,indx_plus_strike) = 1;
                       col(2,indx_plus_dip) = Fdip;
                       smooth = [smooth;col];
%                        disp([indx_minus_strike,indx_plus_strike,indx_minus_dip,indx_plus_dip]);
                   end
               end
           end
       end
   end
end