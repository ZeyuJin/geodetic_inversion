function azi=get_insar_azimuth(X,Y,z);
% 
%  get the azimuth of the flying direction
%  
%  Usage: azi=get_insar_azimuth(X,Y,z);
%
%  X,Y should be geographic coordinates
%
format long
d2r=pi/180;

x=X(1,:);
y=Y(:,1);

h=figure('units','normalized','outerposition',[0 0 1 1]);
h=imagesc(x,y,z);
set(gca,'YDir','Normal');
axis equal
colorbar
colormap('jet');
hold on

Nclick=20;

lon_pt=zeros(Nclick,1);
lat_pt=zeros(Nclick,1);

for k=1:Nclick;
    [lon_pt(k),lat_pt(k)]=ginput(1);
    plot(lon_pt(k),lat_pt(k),'ro','LineWidth',2,'MarkerFaceColor','g','MarkerSize',8);
    title(['# of pts selected: ',num2str(k),'/',num2str(Nclick)]);
end

[lat_sort,indx_sort]=sort(lat_pt);
lon_sort=lon_pt(indx_sort);

lon0=mean(lon_sort);
lat0=mean(lat_sort);
lon_center=lon0;

[xo,yo]=ll2xy(lon0,lat0,lon_center);
[xutm,yutm]=ll2xy(lon_sort,lat_sort,lon_center);
xpt=xutm-xo;
ypt=yutm-yo;


k=polyfit(xpt,ypt,1);
yp=polyval(k,xpt);

xutm_plot=xpt+xo;
yutm_plot=yp+yo;
[lon_plot,lat_plot]=xy2ll(xutm_plot,yutm_plot,lon_center);

plot(lon_plot,lat_plot,'k-','LineWidth',2);

theta=atan(k(1));
alpha=theta/d2r;
if (alpha>0);
 azi=180+90.0-alpha; %descending
else
 azi=-(90.0+alpha); %ascending
end

azi_display=sprintf('%0.2f',azi);
title(['Satellite Heading Azimuth: ',azi_display])
