function [h,h1,h2]=plot_insar_sample(xinsar,yinsar,zinsar,zout,xx1,xx2,yy1,yy2);
% plot out InSAR image and its downsampling
% 
% Usage: h=plot_insar_sample(xinsar,yinsar,zinsar,zout,xx1,xx2,yy1,yy2);
%
%
% by Kang Wang on 08/27/2015

format long
set(0,'defaultTextFontName', 'Hevitica');
set(0,'defaultTextFontSize', 12);

Npt=length(zout);


XX=[xx1';xx2';xx2';xx1'];
YY=[yy1';yy1';yy2';yy2'];
C=[zout';zout';zout';zout'];

zmin=min(min(zinsar(~isnan(zinsar)))); 
zmax=max(max(zinsar(~isnan(zinsar)))); 
rz=zmax-zmin;

cmin=zmin+0.05*rz;
cmax=zmax-0.05*rz;

xmin=min(xinsar);
xmax=max(xinsar);
ymin=min(yinsar);
ymax=max(yinsar);

h=figure('units','normalized','outerposition',[0 0 1 1]);
set(h,'renderer','painters');

psv1=[0.05 0.3 0.45 0.55];
h1=subplot('position',psv1);
h=imagesc(xinsar,yinsar,zinsar);
set(gca,'YDir','Normal');
set(h,'alphadata',~isnan(zinsar));
axis equal;
hcolor=colorbar('north','Position',[0.1 0.25,0.3,0.01]);
colormap('jet')
caxis([cmin cmax]);
axis([xmin xmax ymin ymax])

psv2=[0.5 0.3 0.45 0.55];
h2=subplot('position',psv2);
patch(XX,YY,C);
axis equal
hcolor=colorbar('north','Position',[0.55 0.25,0.3,0.01]);
colormap('jet')
caxis([cmin cmax])
axis([xmin xmax ymin ymax])
title(['#pt=',num2str(Npt)]);






