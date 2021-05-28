function [UE,UN,UZ]=calc_green_3d(data_slip_model,xe,yn,varargin)
% Calculate the Green's Function for InSAR observation using Okada's solution
%
% Usage: G=calc_green_insar(data_slip_model,data_insar);
%
%   data_insar=[xe_insar,yn_insar,los_insar,ve_insar,vn_insar,vz_insar];
%
% by Kang Wang on 07/27/2015
% Last Udpated by Kang Wang on 07/27/2015

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
ss = data_slip_model(:,12);
ds = data_slip_model(:,13);

Npatch=length(iall);
Npara=2*Npatch;

% xe_insar=data_insar(:,1);
% yn_insar=data_insar(:,2);
% zinsar=data_insar(:,3);
% ve_insar=data_insar(:,4);
% vn_insar=data_insar(:,5);
% vz_insar=data_insar(:,6);
% Nobs=length(xe_insar);
Nobs = length(xe);
UE = zeros(Nobs,1);
UN = zeros(Nobs,1);
UZ = zeros(Nobs,1);

G=zeros(Nobs,Npara);
nu=0.25;
HF=1;

if ~isempty(varargin)
    for CC = 1:floor(length(varargin)/2)
        try
            switch lower(varargin{CC*2-1})
                case 'poisson'
                    nu = varargin{CC*2};
            end
        catch
            error('Unrecognized Keyword');
        end
    end
end

for k=1:Npatch
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

%   xpt=xe_insar-xxo;
%   ypt=yn_insar-yyo;
  xpt = xe - xxo;
  ypt = yn - yyo;
  delta=dip0(k)*d2r;
  d=-zzo;
  len=lp(k);
  W=wp(k);
  fault_type1=1;
  fault_type2=2; 
  
%   U1=1;
%   U2=1;
  U1 = ss(k);
  U2 = ds(k);
  tp=zeros(size(xe));
  [ue1,un1,uz1]=calc_okada(HF,U1,xpt,ypt,nu,delta,d,len,W,fault_type1,strike,tp);
  [ue2,un2,uz2]=calc_okada(HF,U2,xpt,ypt,nu,delta,d,len,W,fault_type2,strike,tp);
%   ulos1=ue1.*ve_insar+un1.*vn_insar+uz1.*vz_insar;
%   ulos2=ue2.*ve_insar+un2.*vn_insar+uz2.*vz_insar;
  UE = UE + ue1 + ue2;
  UN = UN + un1 + un2;
  UZ = UZ + uz1 + uz2;

%   G(:,k)=ulos1;
%   G(:,k+Npatch)=ulos2;

end

end