function G=calc_green_AZO_okada(data_slip_model,data_insar);
% Calculate the Green's Function for Azimuth offsets observation using Okada's solution
%
% Usage: G=calc_green_AZO_okada(data_slip_model,data_insar);
%
%   data_insar=[xe_azo,yn_azo,azo,ve_insar,vn_insar,vz_insar];
%
% by Kang Wang on 07/27/2015
% Last Udpated by Kang Wang on 07/27/2015
% Modified by Zeyu Jin on 10/07/2019

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

xe_insar=data_insar(:,1);
yn_insar=data_insar(:,2);
zinsar=data_insar(:,3);
ve_insar=data_insar(:,4);
vn_insar=data_insar(:,5);
vz_insar=data_insar(:,6);
Nobs=length(xe_insar);

% convert the LOS angle to heading angle
theta_az = -atan2d(vn_insar,ve_insar) - 180;

G=zeros(Nobs,Npara);
nu=0.25;
HF=1;

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

  xpt=xe_insar-xxo;
  ypt=yn_insar-yyo;
  delta=dip0(k)*d2r;
  d=-zzo;
  len=lp(k);
  W=wp(k);
  fault_type1=1;
  fault_type2=2; 
  
  U1=1;
  U2=1;
  tp=zeros(size(xe_insar));
  [ue1,un1,uz1]=calc_okada(HF,U1,xpt,ypt,nu,delta,d,len,W,fault_type1,strike,tp);
  [ue2,un2,uz2]=calc_okada(HF,U2,xpt,ypt,nu,delta,d,len,W,fault_type2,strike,tp);
  ulos1=ue1.*sind(theta_az) + un1.*cosd(theta_az);
  ulos2=ue2.*sind(theta_az) + un2.*cosd(theta_az);

  G(:,k)=ulos1;
  G(:,k+Npatch)=ulos2;

end

end
