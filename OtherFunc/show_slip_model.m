function show_slip_model(slip_model,varargin)
format long
%set(0,'defaultAxesFontName', 'AvantGarde')
%set(0,'defaultAxesFontSize', 15)
d2r=pi/180;

% data=load(slip_model);
data = slip_model;
% iflt=data(:,1);
% indx=data(:,2);
% indx_layer=data(:,3);
xe=data(:,4);
yn=data(:,5);
zr=data(:,6);
lp=data(:,7);
wp=data(:,8);
p_strk=data(:,9);
p_dip=data(:,10)*d2r;
tp=data(:,11);
slip1=data(:,12)/100;
slip2=data(:,13)/100;
% min_s1=min(slip1);
% max_s1=max(slip1);
% min_s2=min(slip2);
% max_s2=max(slip2);

% nflt=length(unique(iflt));


% W=sum(unique(wp));
% Ny=length(unique(indx_layer));
N=length(xe);

cmax = 500;
seis_matfile = [];
axis_range = [-30 20 0 50 -25 0];
fault_file = [];
lon_eq = -117.5; 
lat_eq = 35.5;
ref_lon = lon_eq;
% title = 'cm';

if ~isempty(varargin)
   for CC = 1:floor(length(varargin)/2)
       try
          switch lower(varargin{CC*2-1})
              case 'misfit_range'
                  cmax = varargin{CC*2};
              case 'seismic'
                  seis_matfile = varargin{CC*2};
              case 'axis_range'
                  axis_range = varargin{CC*2};
              case 'fault_file'
                  fault_file = varargin{CC*2};
              case 'ref_lon'
                  ref_lon = varargin{CC*2};
              case 'lonc'
                  lon_eq = varargin{CC*2}; 
              case 'latc'
                  lat_eq = varargin{CC*2};
%               case 'title'
%                   title = varargin{CC*2};
          end
       catch
          error('Unrecognized Keyword');
       end
   end
end

figure; hold on
SLIP=[];
XV=[];
YV=[];
ZV=[];
XO=[];
YO=[];
ZO=[];

for i=1:N
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
    
    
    X=[x1,x2,x3,x4]/1000;
    Y=[y1,y2,y3,y4]/1000;
    Z=[z1f,z2f,z3f,z4f]/1000;
    
    xo=mean(X);
    yo=mean(Y);
    zo=mean(Z);
    
    
%     C1=[slip1(i),slip1(i),slip1(i),slip1(i)];
%     C2=[slip2(i),slip2(i),slip2(i),slip2(i)];
    C=sqrt(slip1(i)^2+slip2(i)^2)*ones(1,4);
    patch(X,Y,Z,C,'FaceAlpha',1);
    hold on
    xuf=slip1(i);
    yvf=slip2(i)*cos(p_dip(i));
%     zwf=slip2(i)*sin(p_dip(i))+(slip2(i)*sin(p_dip(i)))/20;
%     maybe a possible bug here
    zwf = slip2(i)*sin(p_dip(i));  
    
    
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
 quiver3(XO,YO,ZO,XV/slipmax,YV/slipmax,ZV/slipmax,1.5,'LineWidth',0.5,'Color',[0,0,0]);

 %set(gca,'HeadStyle','vback1')
  % quiver3d(XO,YO,ZO,XV/slipmax,YV/slipmax,ZV/slipmax,[0,0,0]);
 % lighting phong;
%  camlight head;

if ~isempty(seis_matfile)
    dseis = load(seis_matfile);
    slon = dseis(:,1);
    slat = dseis(:,2);
    sdepth = dseis(:,3) * -1;
    
%     [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);
%     [xx,yy] = utm2ll(slon,slat,0,1);
    [xo,yo] = ll2xy(lon_eq,lat_eq,ref_lon);
    [xx,yy] = ll2xy(slon,slat,ref_lon);
    xs = (xx - xo) ./ 1000;
    ys = (yy - yo) ./ 1000;
    
    xv = [-3 -5 20 20 -2];
    yv = [33 51 51 30 32];
    out = ~inpolygon(xs,ys,xv,yv);
    xin = xs(out);
    yin = ys(out);
    din = sdepth(out);
    
    S = 5;
    scatter3(xin,yin,din,S,'black','filled');   
end

if ~isempty(fault_file)
   for ii = 1:LS
       slon = [lonf(ii) lonf(ii+LS)];
       slat = [latf(ii) latf(ii+LS)];
%        [xx,yy] = utm2ll(slon,slat,0,1);
       [xx,yy] = ll2xy(slon,slat,ref_lon);
       xs = (xx - xo) ./ 1000;
       ys = (yy - yo) ./ 1000;
       line(xs,ys,'color','black','linewidth',1.5);
   end
end

% xlabel('Easting (km)');
% ylabel('Northing (km)');
zlabel('Depth(km)');
axis equal
% zlim([-25 0]);
% axis([-60 -10 285 305 -25 0]);
axis(axis_range);
colormap jet
% hc = colorbar('southoutside');
% hc = colorbar('southoutside');
hc = colorbar;
title(hc,'slip (m)');
if cmax < 5
    caxis([-cmax cmax]);
else
    caxis([0 cmax/100]);
end
% colorbar off
% oldcmap = colormap;
% colormap( flipud(oldcmap) );

grid on
set(gca,'Fontsize',20,'fontweight','bold');
set(gcf,'PaperPositionMode','auto');

end