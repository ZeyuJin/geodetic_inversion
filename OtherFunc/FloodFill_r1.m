%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FloodFill.m implements 2D Goldstein branch cut phase unwrapping algorithm.
% function [IM_unwrapped, rowref, colref] = FloodFill(IM_phase_rads, IM_mag, branch_cuts, IM_mask, colref, rowref)
% It unwraps the phase image, avoiding all branch cuts.
% It can also be used to unwrap phases in cases when there are no branch cuts.
%
% Inputs:
%  IM_phase_rads = 2D array of wrapped phases (rads)
%  IM_mag        = 2D array of magnitudes
%  branch_cuts   = 2D array of identifying branch_cut locations,
%                (Compute using BranchCuts.m or BranchCuts2.m)
% IM_mask is not used here, but it is used in GoldsteinUnwrap2D2.m for the plots.
% colref, rowref = the coordinates of the refernce point,
%                  which is the starting point for the unwrapping algorithm.
% Outputs:
%  IM_unwrapped  = the 2D array of unwrapped phases (rads)
%
% References:
% 1. R. M. Goldstein, H. A. Zebken, and C. L. Werner, “Satellite radar interferometry:
%    Two-dimensional phase unwrapping,” Radio Sci., vol. 23, no. 4, pp. 713–720, 1988.
% 2. D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
%    Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.
%
% Created by B.S. Spottiswoode on 12/10/2004
% 2010/07/23  Carey Smith
%             1. Added colref, rowref as input arguments.
%                Moved the logic to chose the reference point from here
%                GoldsteinUnwrap2D.m, in order to be similar to the
%                Quality Guided routines.
%             2. Implemented Itoh's in-lin method to remove 2*pi jumps
%                (rather than calling unwrap).
%             3. When there are 2 or more valid neighbors, use the one with the 
%                largest magnitude for determining the unwrapped phase.
%                (Previous code re-computed the unwrapped phase for each valid 
%                neighbor, over-writing the previous computation each time.)
%             4. Implemented Content Advisor's recommended usage of & vs. &&
%                and | vs. ||
%             5. Initialzed the output unwrapped version of the phase = nan,
%                rather than zero.  This allows us to distinguish between valid 
%                zeroes and pixels that were not computed.
%                (In plots, the nan's show up as the minimum value.)
%                If you prefer zeros, the original line can be uncommented 
%                or you can reset the nan vaules at the end.
%             6. Allow the edge pixels to be unwrapped.
%             7. Modified the logic so that it shouldn't get "stuck".
%                Check if r_adjoin is empty.  If so, then we are done for that
%                set of points, so break the inner loop.
%                Also, verify that the same set of points are being re-tried,
%                not just that the length of the vector is the same.
%                if there are no valid adjoining pixels, set 
%                adjoin(r_active,c_active) = 0, so that this point isn't re-tried.
%             8. Added logic to skip the branch-cut block, when there are no
%                branch cuts.
%             9. Added additional comments and error-catching logic.
%             10.Changed i to ii and j to jj to avoid redefining the imaginary
%                units
% 2010/07/26  Carey Smith: Made a speed improvement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IM_unwrapped = FloodFill2(IM_phase, IM_mag, branch_cuts, IM_mask, colref, rowref)
% IM_mask is not used

if (branch_cuts(rowref,colref)==1)
  error(['Selected point: (',int2str(rowref),',',int2str(colref),' corresponds to a branch cut.']);
end
n_branch_cuts = length(find(branch_cuts) == 1);
%disp(['FloodFill, number of branch_cuts=',int2str(n_branch_cuts)])

[r_dim, c_dim]=size(IM_phase);
%IM_unwrapped  =zeros(r_dim,c_dim);        % Initialze the output unwrapped version of the phase
IM_unwrapped  = nan(r_dim,c_dim);         % Initialze the output unwrapped version of the phase
% Initialize IM_unwrapped = NaN, so that we can tell which values have been set.
unwrapped_binary=zeros(r_dim,c_dim);      % 0: not yet unwrapped; 1: unwrapped
adjoin =zeros(r_dim,c_dim);               % 0: unwrapped, or not ready; 1: adjoins an unwrapped pixel

% Identify the four pixels adjoining the reference
if(rowref > 1)
  adjoin(rowref-1, colref)=1;
end
if(rowref < r_dim)
  adjoin(rowref+1, colref)=1;
end
if(colref > 1)
  adjoin(rowref, colref-1)=1;
end
if(colref < c_dim)
  adjoin(rowref, colref+1)=1;
end
IM_unwrapped(rowref, colref) = IM_phase(rowref, colref); %Use the reference pixel as the anchor point
unwrapped_binary(rowref, colref) = 1;                    %Mark the first pixel as unwrapped

%disp('Performing floodfill operation ...');
count_limit=0;
r_adjoin = []; % initialize
c_adjoin = []; % initialize
while sum(sum(adjoin))~=0                               %Loop until there are no adjoining pixels or they all lie on the border
  while count_limit<2                                   %or the code gets stuck because of isolated regions
    r_adjoin_stuck = r_adjoin;  % set to previous value
    c_adjoin_stuck = c_adjoin;  % set to previous value
    [r_adjoin, c_adjoin]=find(adjoin);                  %Obtain coordinates of adjoining unwrapped phase pixels
    if(isempty(r_adjoin))
      break;
    end
    % Can get stuck if only edge pixels are adjoining
    % These were avoided in the "if( (r_active<=r_dim-1)...",
    % but were set as adjoining in the loop logic
    if size(r_adjoin)==size(r_adjoin_stuck)
      if(all(r_adjoin==r_adjoin_stuck) && all(c_adjoin==c_adjoin_stuck)) % same points?
        count_limit=count_limit+1;                      %Make sure loop doesn't get stuck
        if(count_limit > 1)
          disp(['FloodFill Caution 1: r_active=',int2str(r_active),', c_active=',int2str(c_active),', count_limit=',num2str(count_limit),' > 1'])
        end
      end
    else
      count_limit=0;
    end
    
    temp=size(r_adjoin);
    for ii=1:temp(1)
      r_active = r_adjoin(ii); % Even tho adjoin gets updated in this loop, r_adjoin & c_adjoin don't
      c_active = c_adjoin(ii);
      phasev   = nan(1,4);     % Initialize.  Will overwrite for valid pixels
      IM_magv  = nan(1,4);     % Initialize.  Will overwrite for valid pixels
        %First search below for an adjoining unwrapped phase pixel
      if(r_active+1<=r_dim)  % Is this a valid index?
        if branch_cuts(r_active+1, c_active)==0
          if unwrapped_binary(r_active+1, c_active)==1
            phase_ref = IM_unwrapped(r_active+1, c_active);       %Obtain the reference unwrapped phase
            % Itoh's Method (suggested by 'Eric' on MATLAB Central to use for a length 2 vector):
            D = IM_phase(r_active, c_active)-phase_ref;
            deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
            phasev(1) = phase_ref + deltap;  % This is the unwrapped phase
            IM_magv(1)= IM_mag(r_active+1, c_active);
          else % unwrapped_binary(r_active+1, c_active)==0
            adjoin(r_active+1, c_active) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
          end
        end
      end
        %Then search above
      if(r_active-1>=1)  % Is this a valid index?
        if branch_cuts(r_active-1, c_active)==0
          if unwrapped_binary(r_active-1, c_active)==1
            phase_ref = IM_unwrapped(r_active-1, c_active);                                   %Obtain the reference unwrapped phase
            D = IM_phase(r_active, c_active)-phase_ref;
            deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
            phasev(2) = phase_ref + deltap;  % This is the unwrapped phase
            IM_magv(2)= IM_mag(r_active-1, c_active);
          else % unwrapped_binary(r_active-1, c_active)==0 
            adjoin(r_active-1, c_active) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
          end
        end
      end
        %Then search on the right
      if(c_active+1<=c_dim)  % Is this a valid index?
        if branch_cuts(r_active, c_active+1)==0
          if unwrapped_binary(r_active, c_active+1)==1
            phase_ref = IM_unwrapped(r_active, c_active+1);                                   %Obtain the reference unwrapped phase
            D = IM_phase(r_active, c_active)-phase_ref;
            deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
            phasev(3) = phase_ref + deltap;  % This is the unwrapped phase
            IM_magv(3)= IM_mag(r_active, c_active+1);
          else % unwrapped_binary(r_active, c_active+1)==0
            adjoin(r_active, c_active+1) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
          end
        end
      end
        %Finally search on the left
      if(c_active-1>=1)  % Is this a valid index?
        if branch_cuts(r_active, c_active-1)==0 
          if unwrapped_binary(r_active, c_active-1)==1
            phase_ref = IM_unwrapped(r_active, c_active-1);                                   %Obtain the reference unwrapped phase
            D = IM_phase(r_active, c_active)-phase_ref;
            deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
            phasev(4) = phase_ref + deltap;  % This is the unwrapped phase
            IM_magv(4)= IM_mag(r_active, c_active-1);
          else % unwrapped_binary(r_active, c_active-1)==0
            adjoin(r_active, c_active-1) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
          end
        end
      end
      idx_del = ~isnan(phasev);
      if(any(idx_del)) % Any valid adjoining pixels?
          % Use the neighbor with the largest magnitude
        IM_max  = max(IM_magv(idx_del));
        idx_max = find((IM_magv >= 0.99*IM_max) & (idx_del==1));
        IM_unwrapped(r_active, c_active) = phasev(idx_max(1));  % Use the first, if there is a tie
        unwrapped_binary(r_active, c_active) = 1;    % Mark the pixel as unwrapped
        adjoin(r_active, c_active) = 0;              % Remove it from the list of adjoining pixels
      else  % no valid adjoining pixels found
        adjoin(r_active,c_active) = 0;  %Remove the current active pixel from the adjoin list
        continue;
      end
    end % for ii=1:temp(1)
    %figure; imagesc(adjoin),       colormap(gray), colorbar, axis square, axis off, title(['Adjoining pixels, ii=',int2str(ii),', count\_limit=',int2str(count_limit)]);
    %figure; imagesc(IM_unwrapped), colormap(gray), colorbar, axis square, axis off, title(['Pixels unwrapped, ii=',int2str(ii),', count\_limit=',int2str(count_limit)]);
  end % while count_limit<100
end % while sum(sum(adjoin(2:r_dim-1,2:c_dim-1)))~=0

if(0)  % DEBUG
  figure;
  imagesc(unwrapped_binary);
  title(['unwrapped\_binary #1']);
  colormap gray;
  colorbar;
end
if(0)  % DEBUG
  figure;
  imagesc(adjoin);
  title(['adjoin #1']);
  colormap gray;
  colorbar;
end  
if(0)  % DEBUG
  figure;
  imagesc(IM_unwrapped);
  title(['FF IM\_unwrapped #1']);
  colormap gray;
  colorbar;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally, fill in the branch cut pixels that adjoin the unwrapped pixels.
% This can be done because the branch cuts actually lie between the pixels,
% and not on top of them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('Filling in branch cuts that border on unwrapped pixels ...');
adjoin=zeros(r_dim, c_dim);  % Clear previous set--just do those along brnach cuts
%Re-load the adjoining pixel matrix with the branch cut values:
if(n_branch_cuts > 0)
  for ii=2:r_dim-1
    for jj=2:c_dim-1
      if branch_cuts(ii,jj)==1 && ...                         %Identify which branch cut pixel borders an unwrapped pixel
          ( (branch_cuts(ii+1,jj)==0 || branch_cuts(ii-1,jj)==0 || branch_cuts(ii,jj-1)==0 || branch_cuts(ii,jj+1)==0) )
        adjoin(ii,jj)=1;
      end
    end
  end
  
  [r_adjoin, c_adjoin]=find(adjoin);  %Obtain coordinates of adjoining unwrapped phase pixels
  temp=size(r_adjoin);
  for ii=1:temp(1)
    r_active=r_adjoin(ii);
    c_active=c_adjoin(ii);
    phasev   = nan(1,4);     % Initialize.  Will overwrite for valid pixels
    IM_magv  = nan(1,4);     % Initialize.  Will overwrite for valid pixels
    %First search below for an adjoining unwrapped phase pixel
    if(r_active+1<=r_dim)  % Is this a valid index?
      if unwrapped_binary(r_active+1, c_active)==1
        phase_ref = IM_unwrapped(r_active+1, c_active);       %Obtain the reference unwrapped phase
        % Itoh's Method (suggested by 'Eric' on MATLAB Central to use for a length 2 vector):
        D = IM_phase(r_active, c_active)-phase_ref;
        deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
        phasev(1) = phase_ref + deltap;  % This is the unwrapped phase
        IM_magv(1)= IM_mag(r_active+1, c_active);
      else % unwrapped_binary(r_active+1, c_active)==0
        adjoin(r_active+1, c_active) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
      end
    end
    %Then search above
    if(r_active-1>=1)  % Is this a valid index?
      if unwrapped_binary(r_active-1, c_active)==1
        phase_ref = IM_unwrapped(r_active-1, c_active);                                   %Obtain the reference unwrapped phase
        D = IM_phase(r_active, c_active)-phase_ref;
        deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
        phasev(2) = phase_ref + deltap;  % This is the unwrapped phase
        IM_magv(2)= IM_mag(r_active-1, c_active);
      else % unwrapped_binary(r_active-1, c_active)==0
        adjoin(r_active-1, c_active) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
      end
    end
    %Then search on the right
    if(c_active+1<=c_dim)  % Is this a valid index?
      if unwrapped_binary(r_active, c_active+1)==1
        phase_ref = IM_unwrapped(r_active, c_active+1);                                   %Obtain the reference unwrapped phase
        D = IM_phase(r_active, c_active)-phase_ref;
        deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
        phasev(3) = phase_ref + deltap;  % This is the unwrapped phase
        IM_magv(3)= IM_mag(r_active, c_active+1);
      else % unwrapped_binary(r_active, c_active+1)==0
        adjoin(r_active, c_active+1) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
      end
    end
    %Finally search on the left
    if(c_active-1>=1)  % Is this a valid index?
      if unwrapped_binary(r_active, c_active-1)==1
        phase_ref = IM_unwrapped(r_active, c_active-1);                                   %Obtain the reference unwrapped phase
        D = IM_phase(r_active, c_active)-phase_ref;
        deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
        phasev(4) = phase_ref + deltap;  % This is the unwrapped phase
        IM_magv(4)= IM_mag(r_active, c_active-1);
      else % unwrapped_binary(r_active, c_active-1)==0
        adjoin(r_active, c_active-1) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
      end
    end
    idx_del = ~isnan(phasev);
    if(any(idx_del)) % Any valid adjoining pixels?
        % Use the neighbor with the largest magnitude
      IM_max  = max(IM_magv(idx_del));
      idx_max = find((IM_magv >= 0.99*IM_max) & (idx_del==1));
      IM_unwrapped(r_active, c_active) = phasev(idx_max(1));  % Use the first, if there is a tie
      unwrapped_binary(r_active, c_active) = 1;      %Mark the pixel as unwrapped
      adjoin(r_active, c_active) = 0;              %Remove it from the list of adjoining pixels
    else  % no valid adjoining pixels found
      %disp(['FloodFill: no valid adjoining pixels found for r_active=',int2str(r_active),', c_active=',int2str(c_active)])
      adjoin(r_active,c_active) = 0;  %Remove the current active pixel from the adjoin list
      continue;
    end
  end % for ii=1:temp(1)

end % if(n_branch_cuts > 0)
%disp('Done!');
