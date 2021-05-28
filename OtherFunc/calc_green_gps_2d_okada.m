function G=calc_green_gps_2d_okada(data_slip_model,data_gps);
%
% Usage:
%     G=calc_green_gps_2d_okada(data_slip_model,data_gps);
%
%     data_gps=[xe_gps,yn_gps,ue,un];

format long
d2r=pi/180;
coln=data_slip_model(:,1);
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
Npara=2*Npatch;


xe_gps=data_gps(:,1);
yn_gps=data_gps(:,2);

Nstn=length(xe_gps);
Nobs=2*Nstn;  %only use the horizontal components
G=zeros(Nobs,Npara);

HF=1;
nu=0.25;

for k=1:Npatch;
   dxf=lp(k)/2;
   dyf=0;
   dzf=0;
   strike=strkp(k)*d2r;
   theta=(90.0-strkp(k))*d2r;
   [dx,dy]=xy2XY(dxf,dyf,-theta);
   dz=dzf;

   xxo=xp(k)+dx;
   yyo=yp(k)+dy;
   zzo=zp(k)+dz;

   xpt=xe_gps-xxo;
   ypt=yn_gps-yyo;
   delta=dip0(k)*d2r;
   d=-zzo;
   len=lp(k);
   W=wp(k);
   fault_type1=1;
   fault_type2=2;
   U1=1;
   U2=1;
   tp=zeros(size(xe_gps));

   [ue1,un1,uz1]=calc_okada(HF,U1,xpt,ypt,nu,delta,d,len,W,fault_type1,strike,tp);
   [ue2,un2,uz2]=calc_okada(HF,U2,xpt,ypt,nu,delta,d,len,W,fault_type2,strike,tp);

   U1_green_gps=[ue1;un1];
   U2_green_gps=[ue2;un2];
   G(:,k)=U1_green_gps;
   G(:,k+Npatch)=U2_green_gps;
end

