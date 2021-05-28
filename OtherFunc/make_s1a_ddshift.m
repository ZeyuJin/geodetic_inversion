clc
clear

PRMFILE='esd_matlab.PRM';
df=load_PRM(PRMFILE,'FREQUENCY_SEPERATION');
dt=load_PRM(PRMFILE,'AZIMUTH_SAMPLING');
nx=load_PRM(PRMFILE,'NX');
ny=load_PRM(PRMFILE,'NY');

c_tmp=2*pi*df*dt;

fid_x=fopen('dd.x','r');
fid_y=fopen('dd.y','r');

fid_fr=fopen('forward.real','r');
fid_fi=fopen('forward.imag','r');
fid_br=fopen('backward.real','r');
fid_bi=fopen('backward.imag','r');

dd_x=fread(fid_x,[nx,ny],'int');
dd_y=fread(fid_y,[nx,ny],'int');

x=dd_x';
y=dd_y';
xx=x(1,:);
yy=y(:,1);
[X,Y]=meshgrid(xx,yy);

fr=fread(fid_fr,[nx,ny],'single');
fi=fread(fid_fi,[nx,ny],'single');
br=fread(fid_br,[nx,ny],'single');
bi=fread(fid_bi,[nx,ny],'single');
fclose(fid_x);
fclose(fid_y);
fclose(fid_fr);
fclose(fid_fi);
fclose(fid_br);
fclose(fid_bi);

FR=fr';
FI=fi';
BR=br';
BI=bi';

CF=FR+i*FI;
CB=BR+i*BI;

nx_filt1=5;
ny_filt1=5;
c_filt1=ones(nx_filt1,ny_filt1)/(nx_filt1*ny_filt1);  


% multilook the forward and backward interferograms before taking the
% double-differencing
E1=conv2(CF,c_filt1,'same'); %multi-looked complex of the forward interferogram
E2=conv2(CB,c_filt1,'same'); %multi-looked complex of the backward interferogram

c_avg=ones(5,5)/25;
a1=conv2(E1.*conj(E2),c_avg,'same');
amp=a1.*conj(a1);
amp1=conv2(E1.*conj(E1),c_avg,'same');
amp2=conv2(E2.*conj(E2),c_avg,'same');
corr_tmp=amp./(amp1.*amp2); 

nx_filt2=15;
ny_filt2=15;
c_filt2=ones(nx_filt2,ny_filt2)/(nx_filt2*ny_filt2);

corr=conv2(corr_tmp,c_filt2,'same');  %low-pass filtering the correlation file

c_dd=E1.*conj(E2);  %double differencing
c_dd_real=real(c_dd);  %real part of the double differencing
c_dd_imag=imag(c_dd);  %imag part of the double differencing


c_real_filt=conv2(c_dd_real,c_filt2,'same');   %low-pass real and imag part of the double differencing
c_imag_filt=conv2(c_dd_imag,c_filt2,'same');

ddphase=atan2(c_imag_filt,c_real_filt);
ddshift=ddphase/c_tmp;
[nl,nc]=size(ddshift);

% decimation
indx_l_ds=1:2:nl;   
indx_c_ds=1:8:nc;

ddshift=ddshift(indx_l_ds,indx_c_ds);
corr=corr(indx_l_ds,indx_c_ds);
X=X(indx_l_ds,indx_c_ds);
Y=Y(indx_l_ds,indx_c_ds);
corr(isnan(corr))=0;
ddshift(isnan(ddshift))=0;

corr_threshold=0.5;
indx_good=find(corr>=corr_threshold);
indx_bad=find(corr<corr_threshold);
ddshift(indx_bad)=NaN;

xgood=X(indx_good);
ygood=Y(indx_good);
zgood=ddshift(indx_good);
wdata=corr(indx_good);
smoothness=40;


[xa,ya,za]=grdread2('a.grd');
% xout=min(xx):200:max(xx);
% yout=min(yy):50:max(yy);
xout=min(xa):100:max(xa);
yout=min(ya):50:max(ya);
zout=xyz2surface(xgood,ygood,zgood,wdata,xout,yout,smoothness);  %fit the sparse points to a smooth surface
grdwrite2(xout,yout,zout,'ddshift_matlab.grd');
quit;
