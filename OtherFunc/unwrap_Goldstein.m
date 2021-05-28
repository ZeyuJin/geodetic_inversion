clc
clear
format long

prmfile='unwrap_goldstein.PRM';
corr_threshold=load_PRM(prmfile,'CORR_THRESHOLD');
interp=load_PRM(prmfile,'INTERPOLATION');

grd_phase=load_PRM(prmfile,'PHASE_DATA');
grd_corr=load_PRM(prmfile,'CORR_DATA');
unwrap_out=load_PRM(prmfile,'UNWRAP_OUT');

[x,y,z1]=grdread2(grd_phase);
[x,y,z2]=grdread2(grd_corr);

[X,Y]=meshgrid(x,y);
X=double(X);
Y=double(Y);
z2(isnan(z2))=0;  % set NaN of correlation to be 0
z1(z2<corr_threshold)=NaN;  %set low correlation pixels to be NaN

%is_interp='n';
is_interp=load_PRM(prmfile,'INTERPOLATION');
if (strcmp(is_interp,'y'));
    indx_good=~isnan(z1);
    xgood=double(X(indx_good));
    ygood=double(Y(indx_good));
    zgood=double(z1(indx_good));
    disp('interpolating the phase with nearest neighbors...');
    zinterp=griddata(xgood,ygood,zgood,X,Y,'nearest');
    im_phase=zinterp;

else
    z1(isnan(z1))=0; %set the NaN pixels to be zeros;
    im_phase=z1;
end

im_mag=z2;
im_mask=ones(size(im_phase));
im_mask(z2<corr_threshold)=0;
mag_max=max(im_mag(:));
disp(['calculating the phase residues start...']);
residue_charge=PhaseResidues_r1(im_phase,im_mask);
disp(['calculating th phase residues finish ...'])
max_box_radius=8;
disp(['calculating the branch cuts start ...']);
branch_cuts=BranchCuts_r1(residue_charge,max_box_radius,im_mask);

disp(['calculating the branch cuts finish ... '])
im_mask(branch_cuts)=0;

im_mag1=im_mag.*im_mask;
%ref_pick='n';
ref_pick=load_PRM(prmfile,'PICK_REFERENCE');
if(strcmp(ref_pick,'y'));
  im_phase_quality=im_mag1;
  minp = im_phase_quality(2:end-1, 2:end-1); minp = min(minp(:));
  maxp = im_phase_quality(2:end-1, 2:end-1); maxp = max(maxp(:));
  figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), colorbar, axis square, axis off; title('Phase quality map');
  %uiwait(msgbox('Select known true phase reference phase point. Black = high quality phase; white = low quality phase.','Phase reference point','modal'));
  uiwait(msgbox('Select known true phase reference phase point. White = high magnitude; Black = low magnitude.','Phase reference point','modal'));
  [xpoint,ypoint] = ginput(1);        %Select starting point for the guided floodfill algorithm
  colref = round(xpoint);
  rowref = round(ypoint);
  close;    
else
  [r_dim, c_dim]=size(im_phase);
  im_mag1(1,:) = 0;                     %Set magnitude of border pixels to 0, so that they are not used for the reference
  im_mag1(r_dim,:) = 0;
  im_mag1(:,1) = 0;
  im_mag1(:,c_dim) = 0;
  [rowrefn,colrefn] = find(im_mag1 >= 0.99*mag_max);
  rowref = rowrefn(1);                  %Choose the 1st point for a reference (known good value)
  colref = colrefn(1); 
end

disp(['doing the unwrapping start ...']);
if(exist('rowref','var'))
  im_unwrapped = FloodFill_r1(im_phase, im_mag, branch_cuts, im_mask, colref, rowref); % Flood fill phase unwrapping
else
  im_unwrapped = FloodFill_r1(im_phase, im_mag, branch_cuts, im_mask); % Flood fill phase unwrapping
end
disp(['doing the unwrapping finish ...']);

im_unwrapped(isnan(z1))=NaN;
grdwrite2(x,y,im_unwrapped,unwrap_out);

im_mask_new=ones(size(im_unwrapped));
im_mask_new(isnan(im_unwrapped))=0;
grdwrite2(x,y,im_mask_new,'mask_branch_cuts.grd');
quit;
