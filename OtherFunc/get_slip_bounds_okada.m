function slip_bounds=get_slip_bounds_okada(slip_model_in,strk_x);
% get the indxes of patches with pre-defined slip
% % Usage: slip_bounds=get_slip_bounds_okada(slip_model_in,strk_x);
% by Kang Wang on 08/21/2015
% Last Updaetd by Kang Wang on 08/21/2015

format long
d2r=pi/180;

data_slip_model=slip_model_in;
icoln=data_slip_model(:,1);
iall=data_slip_model(:,2);
irow=data_slip_model(:,3);
xp=data_slip_model(:,4);
yp=data_slip_model(:,5);
zp=data_slip_model(:,6);
lp=data_slip_model(:,7);
wp=data_slip_model(:,8);
strkp=data_slip_model(:,9);
dip0=data_slip_model(:,10);
Npatch=length(iall);

theta_x=(90.0-strk_x)*d2r;

xxc=zeros(Npatch,1);
yyc=zeros(Npatch,1);
zzc=zeros(Npatch,1);


for i=1:Npatch;
  dxf=lp(i)/2;
  dyf=0;
  dzf=0;
  theta=(90.0-strkp(i))*d2r;
  [dx,dy]=xy2XY(dxf,dyf,-theta);
  xc=xp(i)+dx;
  yc=yp(i)+dy;
  [xxc(i),yyc(i)]=xy2XY(xc,yc,theta_x);
  zzc(i)=zp(i)-(wp(i)/2)*sin(dip0(i)*d2r);
end

depths=unique(irow);
Ndepth=length(depths);
slip_coff=ones(size(iall));

indx_taper=[];
for k=1:Ndepth;
   indx_this_layer=find(irow==depths(k));
   x_this_layer=xxc(indx_this_layer);
   y_this_layer=yyc(indx_this_layer);
   indx_all_this_layer=iall(indx_this_layer);
   [xmin_this_layer,ixmin]=min(x_this_layer);
   [xmax_this_layer,ixmax]=max(x_this_layer);
   indx_xmin=indx_all_this_layer(ixmin);
   indx_xmax=indx_all_this_layer(ixmax);
   indx_taper=[indx_taper;indx_xmin;indx_xmax];
end

%indx_zmax1=find(irow==1);
%indx_zmax2=find(irow==2);
indx_zmin=find(irow==Ndepth);
indx_taper=[indx_taper;indx_zmin];
%indx_taper=[indx_taper;indx_zmin;indx_zmax1];
indx_taper=unique(indx_taper);
slip_coff(indx_taper)=0.001;

slip_bounds=[slip_coff,slip_coff];

%figure;
%scatter3(xp,yp,zp,20,slip_coff,'filled');
%axis equal
%colorbar

