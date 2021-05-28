% Last Updated by Kang Wang on 10/01/2014
function show_slip_okada(slip_model,alphavalue);
% Usage: show_slip_okada(data,alphavalue);
%
%  Example: show_slip_data(slip_model_in,0.5);

% Last Updated by Kang Wang on 12/15/2014
% Last Updated by Kang Wang on 08/24/2015
% Last Update by Kang Wang on 10/05/2015
format long
set(0,'defaultAxesFontName', 'AvantGarde')
set(0,'defaultAxesFontSize', 15)
d2r=pi/180;

%data=load(slip_model);
data=slip_model;
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
min_s1=min(slip1);
max_s1=max(slip1);

min_s2=min(slip2);
max_s2=max(slip2);

nflt=length(unique(iflt));


W=sum(unique(wp));
Ny=length(unique(indx_layer));

N=length(indx);
% figure
SLIP=[];
XV=[];
YV=[];
ZV=[];
XO=[];
YO=[];
ZO=[];
XX=zeros(4,N);
YY=zeros(4,N);
ZZ=zeros(4,N);
CC=zeros(4,N);
for i=1:N;
    p_theta=(90.0-p_strk(i))*d2r;
    
    [x1f,y1f]=xy2XY(xe(i),yn(i),p_theta);
    %y1f=y1f+zr(i)*cos(p_dip(i));
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
    Z=[z1f,z2f,z3f,z4f];
    
    xo=mean(X);
    yo=mean(Y);
    zo=mean(Z);
    
    
    C1=[slip1(i),slip1(i),slip1(i),slip1(i)];
    C2=[slip2(i),slip2(i),slip2(i),slip2(i)];
    C=sqrt(slip1(i)^2+slip2(i)^2)*ones(1,4);
    XX(1,i)=x1;
    XX(2,i)=x2;
    XX(3,i)=x3;
    XX(4,i)=x4;

    YY(1,i)=y1;
    YY(2,i)=y2;
    YY(3,i)=y3;
    YY(4,i)=y4;
   
    ZZ(1,i)=z1f;
    ZZ(2,i)=z2f;
    ZZ(3,i)=z3f;
    ZZ(4,i)=z4f;
  
    CC(1,i)=C(1);
    CC(2,i)=C(2);
    CC(3,i)=C(3);
    CC(4,i)=C(4);

%    patch(X,Y,Z,C)
%    hold on
    xuf=slip1(i);
    yvf=slip2(i)*cos(p_dip(i));
    zwf=slip2(i)*sin(p_dip(i))+(slip2(i)*sin(p_dip(i)))/20;
    
    
    [xu,yu]=xy2XY(xuf,yvf,-p_theta);
    XV=[XV;xu];
    YV=[YV;yu];
    ZV=[ZV;zwf];
    
    XO=[XO;xo];
    YO=[YO;yo];
    ZO=[ZO;zo];
    
    slip0=sqrt(xu^2+yu^2+zwf^2);
    SLIP=[SLIP;slip0];
%    quiver3(xo,yo,zo,xu,yu,zwf,0.5,'LineWidth',2.0,'Color',[0.0,0.0,0.0]);
%     quiver3D([xo,yo,zo],[xu,yu,zwf],'k');
%     lighting phong;
%     camlight head;
end
slipmax=max(SLIP);
%[XV,YV,ZV]
% quiver3(XO,YO,ZO,XV/slipmax,YV/slipmax,ZV/slipmax,1.5,'LineWidth',0.5,'Color',[0,0,0]);
 %set(gca,'HeadStyle','vback1')
  % quiver3d(XO,YO,ZO,XV/slipmax,YV/slipmax,ZV/slipmax,[0,0,0]);
 % lighting phong;
%  camlight head;

%axis equal

figure
if (sqrt(sum(YO.^2))<10);
  patch(XX/1000,ZZ/1000,CC);
  hold on
  quiver(XO/1000,ZO/1000,XV/slipmax,ZV/slipmax,0.8,'LineWidth',0.5,'Color',[0,0,0])

  xlabel('East (km)');
  ylabel('Depth(km)')
  axis equal
  colorbar
  grid on
elseif (sqrt(sum(XO.^2))<10);
  patch(YY/1000,ZZ/1000,CC);
  hold on
  quiver(YO/1000,ZO/1000,YV/slipmax,ZV/slipmax,0.8,'LineWidth',0.5,'Color',[0,0,0]);
  xlabel('North (km)');
  ylabel('Depth(km)')
  axis equal
  colorbar
  grid on
elseif (sqrt(sum(ZO.^2))<10);
  patch(YY/1000,ZZ/1000,CC);
  hold on
  quiver(XO/1000,YO/1000,XV/slipmax,YV/slipmax,0.8,'LineWidth',0.5,'Color',[0,0,0]);
  xlabel('East (km)');
  ylabel('North (km)')
  axis equal
  colorbar
  grid on

else
  patch(XX/1000,YY/1000,ZZ/1000,CC);
  hold on
  quiver3(XO/1000,YO/1000,YO/1000,XV/slipmax,YV/slipmax,ZV/slipmax,0.8,'LineWidth',0.5,'Color',[0,0,0]);
  xlabel('East (km)');
  ylabel('North (km)');
  ylabel('Depth (km)');
  axis equal
  colorbar
  grid on

end
axis tight



