%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BranchCuts.m generates branch cuts based on the phase residues. This is
% done using the Goldstein method, as described in "Two-dimensional phase
% unwrapping: theory, algorithms and software" by Dennis Ghiglia and
% Mark Pritt.
% "residue_charge" is a matrix wherein positive residues are 1 and
% negative residues are 0.
% "max_box_radius" defines the maximum search radius for the balancing of
% residues. If this is too large, areas will be isolated by the branch
% cuts.
% "IM_mask" is a binary matrix. This serves as an artificial border for the
% branch cuts to connect to.
% Created by B.S. Spottiswoode on 15/10/2004
% Last modified on 18/10/2004
% 07/19/2010 Modified by Carey Smith
%             Corrected a NaN divide by zero problem
%             Return immediately, if no residues
%             Eliminated code Analyzer warnings about && vs. &, || vs. |
% 2010/09/15  Match with residues of opposite sign.
%             Change the order of the loops: Put the radius loop on the outside,
%             to let each residue try to find its closest neighbor of opposite
%             sign.  This generally finds closer pairs.
%             Use residue_charge_masked in place of residue_charge in some places.
%             Modified "for n" loop to just use the needed indices.
%             When making a branch cut to the edge, go straight, rather than diagonally.
%             Eliminated unnecessary sections of code & resulting unused variables.
%             Sped-up internal function "branch_cuts".
%             Added some error checks.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function branch_cuts=BranchCuts(residue_charge, max_box_radius, IM_mask)

branch_cuts = ~IM_mask;                             %Define initial branch cuts borders as the mask.

residue_charge_masked=residue_charge;
residue_charge_masked(logical(~IM_mask))=0;         %Remove all residues except those in the mask
[rowres,colres] = find(residue_charge_masked~=0);             %Find the coordinates of the residues
if(isempty(rowres));
  disp(['BranchCuts: No residues, length(rowres)=',int2str(length(rowres)),...
    '; sum(abs(residue_charge))=',int2str(sum(abs(residue_charge(:)))),...
    '; sum(abs(residue_charge_masked))=',int2str(sum(abs(residue_charge_masked(:))))])
  return;
end                   % no residues

% Allocate spcae
[rowdim, coldim] =size(residue_charge);
residue_balanced =zeros(rowdim, coldim);            %Initially designate all residues as unbalanced

%disp('Calculating branch cuts ...');
%tic;
n_residues = length(rowres);
max_box_radius = min(max_box_radius,floor(length(residue_charge)/2));
for(radius=1:max_box_radius) % Loop thru box sizes
  %disp(['BranchCuts.m: radius=',num2str(radius)])
  for i=1:n_residues;                               %Loop through the residues
    r_active=rowres(i);                             %Coordinates of the active residue
    c_active=colres(i);
    if(residue_balanced(r_active, c_active) > 0)    % Already balanced, possibly by a preceeding residue?
      %disp(['i=',int2str(i),'; r_active=',int2str(r_active),'; c_active=',int2str(c_active),...
      %  '; residue_charge=',int2str(residue_charge(r_active, c_active)),...
      %  '; residue_charge_masked=',int2str(residue_charge_masked(r_active, c_active)),'; Already balanced or not in the mask'])
      continue;
    end
    if r_active<=1 || r_active>=rowdim || c_active<=1 || c_active>=coldim %Is this residue on the image border?
      branch_cuts(r_active, c_active)=1;  % Make this point a branchcut to the edge
      residue_balanced(r_active, c_active) = 1;       %Mark this residue as balanced
      %disp(['i=',int2str(i),'; r_active=',int2str(r_active),'; c_active=',int2str(c_active),'; residue_charge=',int2str(residue_charge(r_active, c_active)),'; edge residue, so balanced'])
      residue_charge_masked(r_active, c_active) = 0;  % Remove from the set of unmatched residues
      continue;
    end
    charge_counter = residue_charge_masked(r_active, c_active);  %Store the initial residue charge
    %disp(['i=',int2str(i),'; r_active=',int2str(r_active),'; c_active=',int2str(c_active),...
    %  '; residue_charge=',int2str(residue_charge(r_active, c_active)),'; radius=',num2str(radius)])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This portion of code serves to search the box perimeter,
    %place branch cuts, and keep track of the summed residue charge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m1 = max(r_active-radius,1);
    m2 = min(r_active+radius,rowdim);
    n1 = max(c_active-radius,1);
    n2 = min(c_active+radius,coldim);
    for m = m1:m2  %Coordinates of the box border pixels
      if(m==m1 || m==m2)  %Ensure that only the border pixels are being tested
        ndel = 1;       %Use all pixels in top & bottom rows
      else              % Avoid re-testing pixels that have already been tested
        ndel = n2-n1;   %Use only the end pixels in middle rows
      end
      for n = n1:ndel:n2  % Coordinates of only the box border pixels
        %disp([' m=',int2str(m),'; n=',int2str(n),...
        %'; residue_charge(m,n)=',int2str(residue_charge(m,n)),...
        %'; residue_charge_masked(m,n)=',int2str(residue_charge_masked(m,n)),...
        %'; IM_mask=',num2str(IM_mask(m,n),2),...
        %'; res_balanced=',int2str(residue_balanced(m,n))])
        if (charge_counter~=0)  % Not yet resolved
          if m<=1 || m>=rowdim          %Is the current pixel on the image border?
            branch_cuts = PlaceBranchCutInternal(branch_cuts, r_active, c_active, m, c_active);  %Place a branch cut between the active point and the mask border
            residue_balanced(r_active, c_active)=1;         %Mark the centre residue as balanced
            %disp([' Balanced to a row edge, branch_cuts(',int2str(r_active),',',int2str(n),')=',int2str(branch_cuts(r_active, n))])
            charge_counter=0;                               %Label the charge as balanced
            residue_charge_masked(r_active, c_active) = 0;  % Remove from the set of unmatched residues
            residue_charge_masked(r_active, n)        = 0;  % Remove from the set of unmatched residues
            break;
          elseif n<=1 || n>=coldim         %Is the current pixel on the image border?
            branch_cuts = PlaceBranchCutInternal(branch_cuts, r_active, c_active, r_active, n);  %Place a branch cut between the active point and the nearest mask border
            residue_balanced(r_active, c_active)=1;         %Mark the centre residue as balanced
            %disp([' Balanced to a col edge'])
            charge_counter=0;                               %Label the charge as balanced
            residue_charge_masked(r_active, c_active) = 0;  % Remove from the set of unmatched residues
            residue_charge_masked(m,        c_active) = 0;  % Remove from the set of unmatched residues
            break;
            % Logic to connect positive residues to negative residues.
          elseif(residue_charge_masked(r_active, c_active) * residue_charge_masked(m,n) == -1) % Is the current pixel a residue of opposite sign?
            branch_cuts = PlaceBranchCutInternal(branch_cuts, r_active, c_active, m, n);  %Place a branch cut regardless of the charge_counter value
            %disp([' m=',int2str(m),'; n=',int2str(n),'; residue_charge_masked(m,n)=',int2str(residue_charge_masked(m,n)),'charge_counter=',int2str(charge_counter),', so balanced'])
            charge_counter=0;                               %Label the charge as balanced
            residue_balanced(r_active, c_active)=1;         %Mark the centre (active) residue as balanced
            residue_balanced(m,        n)       =1;        %Mark the end-pt. residue as balanced
            residue_charge_masked(r_active, c_active) = 0; % Remove from the set of unmatched residues
            residue_charge_masked(m,        n)        = 0; % Remove from the set of unmatched residues
            break;
            %else
            %  disp([' m=',int2str(m),'; n=',int2str(n),'; residue_charge=',int2str(residue_charge(m,n)),...
            %    '; residue_charge_masked(m,n)=',int2str(residue_charge_masked(m,n)),'; charge_counter=',int2str(charge_counter),';not yet balanced'])
          end % if [m<=1 || m>=rowdim || n<=1 || n>=coldim]; elseif(residue_charge_masked(r_active, c_active) * residue_charge_masked(m,n) == -1)
        end % if charge_counter~=0
        if(charge_counter==0); break; end;
      end % for n
      if(charge_counter==0); break; end;
    end % for m
    
  end % for i=1:n_residues
  %disp([' sum(residue_balanced(:))=',int2str(sum(residue_balanced(:)))])
  if(sum(residue_balanced(:)) >= n_residues)  % are all of the residues balanced?
    break;
  end
end % for(radius=1:floor(length(residue_charge)/2))

%t=toc;
%disp(['Branch cut operation completed in ',int2str(t),' seconds.']);
%disp(['Residues: total= ',int2str(n_residues),', accounted= ',int2str(sum(residue_balanced(:))),', unaccounted= ',int2str(n_residues-sum(residue_balanced(:)))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PlaceBranchCutInternal.m places a branch cut between the points [r1, c1] and
% [r2, c2]. The matrix branch_cuts is binary, with 1's depicting a
% branch cut.
%
% 2010/09/15  Carey Smith  Changed atan() to atan2() & reduced to one "for loop"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function branch_cuts = PlaceBranchCutInternal(branch_cuts, r1, c1, r2, c2)

branch_cuts(r1,c1)=1;               %Fill the starting points
branch_cuts(r2,c2)=1;               %Fill the ending points
radius=sqrt((r2-r1)^2 + (c2-c1)^2); %Distance between points
warning off MATLAB:divideByZero;    %This warning has been removed from MATLAB 2010
%theta=atan2((c2-c1),(r2-r1));       % CS: atan2() gives the full line angle (=/-pi),
%                                   %     so, we need just one "for" loop
%cos_theta = cos(theta)             % CS: Compute once, outside the for loop
%sin_theta = sin(theta)             % CS: Compute once, outside the for loop
% Compute cos & sin by trig, w/out computing the angle
dist = sqrt((c2-c1)^2 + (r2-r1)^2);
cos_theta = (r2-r1)/dist;           % CS: Compute once, outside the for loop
sin_theta = (c2-c1)/dist;           % CS: Compute once, outside the for loop
for i=1:radius                      %Number of points to fill in
  r_fill=r1 + round(i*cos_theta);
  c_fill=c1 + round(i*sin_theta);
  
  if(c_fill==0)
    error(['PlaceBranchCutInternal:invalid_c',' c_fill==0',...
      '; r1=',int2str(r1),'; r2=',int2str(r2),'; c1=',int2str(c1),'; c2=',int2str(c2),...
      '; i=',int2str(i),'; cos_theta=',num2str(cos_theta,4),'; sin_theta=',num2str(sin_theta,4)])
  end
  if(r_fill==0)
    error(['PlaceBranchCutInternal:invalid_r',' r_fill==0',...
      '; r1=',int2str(r1),'; r2=',int2str(r2),'; c1=',int2str(c1),'; c2=',int2str(c2),...
      '; i=',int2str(i),'; cos_theta=',num2str(cos_theta,4),'; sin_theta=',num2str(sin_theta,4)])
  end
  
  branch_cuts(r_fill, c_fill) = 1;
end
return;

