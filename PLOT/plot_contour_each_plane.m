function [data,ntour] = plot_contour_each_plane(slip_model,varargin)
% slip_model should correspond to each slip plane
format long
%set(0,'defaultAxesFontName', 'AvantGarde')
%set(0,'defaultAxesFontSize', 15)
d2r=pi/180;

cmax = 10;
axis_range = [0 25 -15 0];
title_name = '';
if ~isempty(varargin)
   for CC = 1:floor(length(varargin)/2)
       try
          switch lower(varargin{CC*2-1})
              case 'misfit_range'
                  cmax = varargin{CC*2};
              case 'axis_range'
                  axis_range = varargin{CC*2};
              case 'title'
                  title_name = varargin{CC*2};
          end
       catch
          error('Unrecognized Keyword');
       end
   end
end

fin = slip_model(:,1);
fault_id = unique(fin);
fault_num = length(fault_id);

for kk = 1:fault_num

this_fault = find(fin == fault_id(kk));
data = slip_model(this_fault,:);
% data = slip_model;
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
% tp=data(:,11);
slip1=data(:,12);
slip2=data(:,13);

% min_s1=min(slip1);
% max_s1=max(slip1);
% min_s2=min(slip2);
% max_s2=max(slip2);
% nflt=length(unique(iflt));

% W=sum(unique(wp));
% Ny=length(unique(indx_layer));
N=length(indx);

figure; hold on
SLIP=[];
XV=[];
YV=[];
ZV=[];
XO=[];
YO=[];
ZO=[];

minlp = min(lp);

for i=1:N
    p_theta=(90.0-p_strk(i)+180)*d2r;
    [xr,~] = xy2XY(xe,yn,p_theta);
    xend = min(xr);
%     disp(xend);
    
    [x1f,y1f]=xy2XY(xe(i),yn(i),p_theta);
%     disp([x1f,y1f]);
    %y1f=y1f+zr(i)*cos(p_dip(i));
    x1f=x1f-xend-(lp(i)-minlp);
    z1f=zr(i) ./ sin(pi - p_dip(i));
%     disp([x1f,y1f]);
    
    x2f=x1f+lp(i);
    y2f=y1f;
    z2f=z1f;
    
    x3f=x2f;
    y3f=y2f;
    z3f=z2f-wp(i);
%     y3f=y2f-wp(i)*cos(p_dip(i));
%     z3f=z2f-wp(i)*sin(p_dip(i));
    
    x4f=x1f;
    y4f=y3f;
    z4f=z3f;
    
%     [x1,y1]=xy2XY(x1f,y1f,-p_theta);
%     [x2,y2]=xy2XY(x2f,y2f,-p_theta);
%     [x3,y3]=xy2XY(x3f,y3f,-p_theta);
%     [x4,y4]=xy2XY(x4f,y4f,-p_theta);
    
    [x1,y1]=xy2XY(x1f,y1f,0);
    [x2,y2]=xy2XY(x2f,y2f,0);
    [x3,y3]=xy2XY(x3f,y3f,0);
    [x4,y4]=xy2XY(x4f,y4f,0);    
    
    X=[x1,x2,x3,x4]/1000;
    Y=[y1,y2,y3,y4]/1000;
    Z=[z1f,z2f,z3f,z4f]/1000;
    
    xo=mean(X);
    yo=mean(Y);
    zo=mean(Z);
    
    
%     C1=[slip1(i),slip1(i),slip1(i),slip1(i)];
%     C2=[slip2(i),slip2(i),slip2(i),slip2(i)];
    C=sqrt(slip1(i)^2+slip2(i)^2)*ones(1,4);
%     patch(X,Y,Z,C,'FaceAlpha',1);
    patch(X,Z,C,'FaceAlpha',1);
    hold on
    xuf=slip1(i);
    yvf=slip2(i)*cos(p_dip(i));
%     zwf=slip2(i)*sin(p_dip(i))+(slip2(i)*sin(p_dip(i)))/20;
    zwf = slip2(i)*sin(p_dip(i));
    
%     [xu,yu]=xy2XY(xuf,yvf,-p_theta);
    [xu,yu]=xy2XY(xuf,yvf,0);
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

% interpolate the whole grid plane
F = scatteredInterpolant(XO,ZO,SLIP);
xq = linspace(min(XO), max(XO), ceil(length(indx)/max(indx_layer))*2);
yq = linspace(min(ZO), max(ZO), max(indx_layer)*2);
[xq,yq] = meshgrid(xq,yq);
VQ = F(xq,yq);

% add contour plot
% slipmax=max(SLIP);
ntour = floor(max(max(VQ))/50);   % find the previous cmax bug here
contour(xq,yq,VQ,ntour,'w','LineWidth',3);  % convert to meters

data = [xq,yq,VQ];

% % quiver plot
% slipmax=max(SLIP);
% quiver(XO,ZO,slip1./slipmax,slip2./slipmax,1.5,'LineWidth',0.5,'Color',[0,0,0]);

%  quiver3(XO,YO,ZO,XV/slipmax,YV/slipmax,ZV/slipmax,1.5,'LineWidth',0.5,'Color',[0,0,0]);
 %set(gca,'HeadStyle','vback1')
  % quiver3d(XO,YO,ZO,XV/slipmax,YV/slipmax,ZV/slipmax,[0,0,0]);
 % lighting phong;
%  camlight head;

% xlabel('Easting (km)');
% ylabel('Northing (km)');
% zlabel('Depth(km)')
axis equal
% axis([-30 20 0 50 -25 0]);
% xlim([0 25]);
% ylim([-15 0]);
% title('Coseismic Slip on the NW branch');
title(title_name);
axis(axis_range);
colormap jet
colorbar
caxis([0 cmax]);
grid on
set(gca,'Fontsize',25);
set(gcf,'PaperPositionMode','auto');

end

end