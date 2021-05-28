%filtering the double-differencing interferograms
%by Kang Wang in April. 2016
clc
clear

PRMFILE='dd_matlab.PRM';

Nx=load_PRM(PRMFILE,'NX');
Ny=load_PRM(PRMFILE,'NY');

nx=Nx;
ny=floor(Ny/2);

%Ny=floor(ny/2);
%Ny=ny;
%dir_int='./azi_int/';

dir_int='./';
fid_fr=fopen([dir_int,'forward.real'],'r');
fid_fi=fopen([dir_int,'forward.imag'],'r');
fid_br=fopen([dir_int,'backward.real'],'r');
fid_bi=fopen([dir_int,'backward.imag'],'r');

fr=fread(fid_fr,[ny,nx],'single');
fi=fread(fid_fi,[ny,nx],'single');
br=fread(fid_br,[ny,nx],'single');
bi=fread(fid_bi,[ny,nx],'single');

fclose(fid_fr);
fclose(fid_fi);
fclose(fid_br);
fclose(fid_bi);

% FR=fr(1:Ny,:);
% FI=fi(1:Ny,:);
% BR=br(1:Ny,:);
% BI=bi(1:Ny,:);

FR=fr;
FI=fi;
BR=br;
BI=bi;

% nx_look=10;
% ny_look=10;
% 
% c_filt1=ones(nx_look,ny_look)/(nx_look*ny_look);
% FR_filt=conv2(FR,c_filt1,'same');
% FI_filt=conv2(FI,c_filt1,'same');
% 
% BR_filt=conv2(BR,c_filt1,'same');
% BI_filt=conv2(BI,c_filt1,'same');

FR_filt=imgaussfilt(FR,[4,4]);
FI_filt=imgaussfilt(FI,[4,4]);

BR_filt=imgaussfilt(BR,[4,4]);
BI_filt=imgaussfilt(BI,[4,4]);

clear FR FI BR BI

E1=FR_filt+i*FI_filt;
E2=BR_filt+i*BI_filt;

clear FR_filt FI_filt BR_filt BI_filt

c_dd=E1.*conj(E2);  %double differencing
c_dd_real=real(c_dd);  %real part of the double differencing
c_dd_imag=imag(c_dd);  %imag part of the double differencing

clear c_dd 
dd_real_filt=imgaussfilt(c_dd_real,[5,5]);
dd_imag_filt=imgaussfilt(c_dd_imag,[5,5]);

clear c_dd_real c_dd_image

ddphase=atan2(dd_imag_filt,dd_real_filt);
nx_look=4;
ny_look=4;
filt_look=ones(4,4)/(nx_look*ny_look);
ddphase_mlook=conv2(ddphase,filt_look,'same');
[nl,nc]=size(ddphase);
xin=linspace(0,Nx,nc);
yin=linspace(0,Ny,nl);

indx_y_out=1:ny_look:nl;
indx_x_out=1:nx_look:nc;

xout=xin(indx_x_out);
yout=yin(indx_y_out);
zout=ddphase_mlook(indx_y_out,indx_x_out);
clear ddphase ddphase_mlook
grdwrite2(xout,yout,zout,'ddphase_ra.grd');
quit;


 
