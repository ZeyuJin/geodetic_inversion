%interpolate the phase with nearest neighbor for snaphu unrapping
% by Kang Wang on Mar. 7th. 2016

clc
clear
format long

PRMFILE='snaphu_matlab.PRM';
snaphu_threshold=load_PRM(PRMFILE,'SNAPHU_THRESHOLD');
phase_in = load_PRM(PRMFILE,'PHASE_IN');
corr_in =load_PRM(PRMFILE,'CORR_IN');

phase_out=load_PRM(PRMFILE,'PHASE_OUT');
corr_out=load_PRM(PRMFILE,'CORR_OUT');

%mask_out=load_PRM(PRMFILE,'MASK_OUT');
%snaphu_threshold=0.12;
%[x1,y1,z1]=grdread2('phase_merge_ll.grd');
%[x2,y2,z2]=grdread2('corr_merge_ll.grd');
[x1,y1,z1]=grdread2(phase_in);
[x2,y2,z2]=grdread2(corr_in);

[X,Y]=meshgrid(x1,y1);

z1=double(z1);
X=double(X);
Y=double(Y);

indx_bad_corr=find(z2<snaphu_threshold); %get the index of bad correlation pixels
indx_good_corr=find(z2>=snaphu_threshold);

z1(indx_bad_corr)=NaN;
z1(z1>pi|z1<-pi)=NaN;
indx_good=~isnan(z1);

xgood=X(indx_good);
ygood=Y(indx_good);
zgood=z1(indx_good);

z1_interp=griddata(xgood,ygood,zgood,X,Y,'nearest');
z1_interp(isnan(z2))=0;

z2_interp=z2;
z2_interp(indx_bad_corr)=0;
z2_interp(isnan(z2))=0;

zmask=NaN(size(z2));
zmask(indx_good_corr)=1;

%grdwrite2(x1,y1,z1_interp,'phase_merge_ll_interp.grd');
%grdwrite2(x2,y2,z2_interp,'corr_merge_ll_interp.grd');
%grdwrite2(x2,y2,zmask,'mask_ll_unwrap.grd');

grdwrite2(x1,y1,z1_interp,phase_out);
grdwrite2(x2,y2,z2_interp,corr_out);
%grdwrite2(x2,y2,zmask,mask_out);

quit;
