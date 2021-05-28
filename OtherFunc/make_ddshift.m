% by Kang Wang on 08/29/2016
% Last Updated by Kang Wang on 08/30/2016
% Last Updated by Kang Wang on 09/01/2016
%
%

clc
clear

set(0,'defaultTextFontName', 'Helvetica')
set(0,'defaultTextFontSize',12)

prmfile='ddshift.PRM';
data=load('ddphase');

x=data(:,1);
y=data(:,2);
zphase=data(:,3);
corr=data(:,4);

df = load_gmtsar_PRM(prmfile,'Frequency_separation'); 
dt = load_gmtsar_PRM(prmfile,'Azimuth_sampling');
burst_offset=load_gmtsar_PRM(prmfile,'Burst_offset');
%y=y+burst_offset; %no need any more since the output ddshift has now included this shift

z=zphase/(2*pi*df*dt);
xmin=0;
xmax=load_gmtsar_PRM(prmfile,'num_rng_bins');

ymin=0;
ymax=load_gmtsar_PRM(prmfile,'num_lines');

xout=xmin:100:xmax;
yout=ymin:25:ymax;

[XOUT,YOUT]=meshgrid(xout,yout);


%wdata=corr;           %weighting using correlation
wdata=ones(size(z));   %uniform weighting
smoothness=20;
zout=xyz2surface(x,y,z,wdata,xout,yout,smoothness);

hf1=figure;
set(hf1,'visible','off');
imagesc(xout,yout,zout);
set(gca,'YDir','Normal')
colorbar
colormap('jet');
caxis([-0.02 0.02])
saveas(hf1,'ddshift','eps2c');
%
%

grdwrite2(xout,yout,zout,'ddshift_matlab.grd');
quit;
