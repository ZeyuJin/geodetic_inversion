function [zres,zramp]=get_ifg_ramp(X,Y,Z,ndeg,N);
%
% estimate the linear ramp in an interferogram
%
% Usage: [zres,zramp]=get_ifg_ramp(X,Y,Z,ndeg,N);
%
% X,Y,Z ----- matrixes of the same size
%
% ndeg --- maximum degrees of polynomials
%
% N  --- number of polygons to estimate the ramp
%
% by Kang Wang in Nov. 2015


%x=X(1,:);
%y=Y(:,1);
%z=Z;

%lon0=mean(x);
%lat0=mean(y);
%lon_center=lon0;

%[xo,yo]=ll2xy(lon0,lat0,lon_center);
%[xutm1,yutm1]=ll2xy(x',y(1)*ones(size(x')),lon_center);
%xe1=xutm1-xo;
%yn1=yutm1-yo;

%[xutm2,yutm2]=ll2xy(x(1)*ones(size(y)),y,lon_center);
%xe2=xutm2-xo;
%yn2=yutm2-yo;

%[XX,YY]=meshgrid(xe1,yn2);

%xx=xe1;
%yy=yn2;

XX=X;
YY=Y;
z=Z;
xx=XX(1,:);
yy=YY(:,1);
xmin=min(xx);
xmax=max(xx);
ymin=min(yy);
ymax=max(yy);

indx_good=~isnan(z);
xgood=XX(indx_good);
ygood=YY(indx_good);
zgood=z(indx_good);


figure('units','normalized','outerposition',[0 0 1 1]);
h=imagesc(xx,yy,z);
set(gca,'YDir','Normal');
axis equal
set(h,'alphadata',~isnan(z));
colorbar
colormap('jet')
axis([xmin xmax ymin ymax]);
hold on
Nclick=6;

XFIND=[];
YFIND=[];
ZFIND=[];
for n=1:N;
  lon_pt=zeros(Nclick,1);
  lat_pt=zeros(Nclick,1);
  for k=1:Nclick;
    [lon_pt(k),lat_pt(k)]=ginput(1);
    plot(lon_pt(k),lat_pt(k),'ro','LineWidth',2,'MarkerFaceColor','g','MarkerSize',8);
  end
  plot([lon_pt;lon_pt(1)],[lat_pt;lat_pt(1)],'k-','LineWidth',2);
  in_side=inpolygon(xgood,ygood,lon_pt,lat_pt);
  xfind=xgood(in_side);
  yfind=ygood(in_side);
  zfind=zgood(in_side);
  XFIND=[XFIND;xfind];
  YFIND=[YFIND;yfind];
  ZFIND=[ZFIND;zfind];
end


model=polyfitn([xfind,yfind],zfind,ndeg);
zramp=polyvaln(model,[XX(:),YY(:)]);
zramp=reshape(zramp,size(XX));
zres=z-zramp;
%[a1,a2,a3]=xyz_trend(xfind,yfind,zfind);
%zramp=XX*a1+YY*a2+a3;
%zres=z-zramp;

%[a1,a2,a3,a4]=xyz_trend_bilinear(xfind,yfind,zfind);
%zramp=a1*XX.*YY+a2*XX+a3*YY+a4;
%zres=z-zramp;

function [a1,a2,a3]=xyz_trend(x,y,z);
% Estimate a planar ramp from give data points
% Usage: [a1,a2,a3]=xyz_trend(x,y,z)
%

A=[x,y,ones(size(x))];
b=pinv(A)*z;

a1=b(1)
a2=b(2)
a3=b(3)
           

function [a1,a2,a3,a4,a5,a6]=xyz_trend_bilinear(x,y,z);
%A=[x,y,x.*y,ones(size(x))];
A=[x.*y,x,y,ones(size(x))];
b=pinv(A)*z;

a1=b(1)
a2=b(2)
a3=b(3)
a4=b(4)

