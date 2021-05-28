function [xout,yout,zout]=slip2grid_ll(slip_model,lon0,lat0,lon_center);
% 
% Usage: [xout,yout,zout]=slip2grid_ll(slip_model,lon0,lat0,lon_center);
%
% Last Updated by Kang Wang on July 14, 2016

d2r=pi/180;
% data=slip_model;
data=load(slip_model);
iflt=data(:,1);
indx=data(:,2);
indx_layer=data(:,3);
xe=data(:,4);
yn=data(:,5);
zr=data(:,6);
lp=data(:,7);
wp=data(:,8);
p_strk=data(:,9);
p_dip=data(:,10)*d2r;
tp=data(:,11);
slip1=data(:,12);
slip2=data(:,13);

N=length(indx);
xin=[];
yin=[];
zin=[];

for i=1:N;
    p_theta=(90.0-p_strk(i))*d2r;
    [x1f,y1f]=xy2XY(xe(i),yn(i),p_theta);
   z1f=zr(i);

    x2f=x1f+lp(i);
    y2f=y1f;
    z2f=z1f;

    x3f=x2f;
    y3f=y2f-wp(i)*cos(p_dip(i));
    z3f=z2f-wp(i)*sin(p_dip(i));

    x4f=x1f;
    y4f=y3f;
    z4f=z3f;

    [x1,y1]=xy2XY(x1f,y1f,-p_theta);
    [x2,y2]=xy2XY(x2f,y2f,-p_theta);
    [x3,y3]=xy2XY(x3f,y3f,-p_theta);
    [x4,y4]=xy2XY(x4f,y4f,-p_theta);

    X=[x1,x2,x3,x4];
    Y=[y1,y2,y3,y4];
    xo=mean(X);
    yo=mean(Y);

    s_patch=sqrt(slip1(i)^2+slip2(i)^2);
    xin=[xin;xo];
    yin=[yin;yo];
    zin=[zin;s_patch];

end

xmin=min(xin);
xmax=max(xin);
ymin=min(yin);
ymax=max(yin);

dx=2e3;
dy=2e3;

xout=xmin:dx:xmax;
yout=ymin:dy:ymax;

[X,Y]=meshgrid(xout,yout);
zout=griddata(xin,yin,zin,X,Y);
[xo,yo]=ll2xy(lon0,lat0,lon_center);

[m,n]=size(X);
Npt=m*n;
xpt=reshape(X,Npt,1);
ypt=reshape(Y,Npt,1);

xpt_utm=xpt+xo;
ypt_utm=ypt+yo;

[lon_pt,lat_pt]=xy2ll(xpt_utm,ypt_utm,lon_center);
LON=reshape(lon_pt,size(X));
LAT=reshape(lat_pt,size(Y));
xout=LON(1,:);
yout=LAT(:,1);

