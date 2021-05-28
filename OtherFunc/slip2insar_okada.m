function zout=slip2insar_okada(XIN,YIN,zin,look_e,look_n,look_z,slip_model_in,varargin)
%
%  Usage: 
%     los_out=slip2insar_okada(XIN,YIN,zin,los_in,look_e,look_n,look_z,slip_model_in);
%    
%%%%% input
%      XIN     matrix of x-coordinates
%      YIN     matrix of y-coordinates
%      zin     matrix of input los values (could cantain NaN)
%      look_e  matrix of eastard unit looking
%      look_n
%      look_z
%
% insar and slip model should in the same cartisian coordinates

% by Kang Wang on 07/29/2015


nu=0.25;
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

format long
d2r=pi/180;
indx_patch=slip_model_in(:,2);
xp=slip_model_in(:,4);
yp=slip_model_in(:,5);
zp=slip_model_in(:,6);
lp=slip_model_in(:,7);
wp=slip_model_in(:,8);
strkp=slip_model_in(:,9);
dip0=slip_model_in(:,10);
s1_patch=slip_model_in(:,12);
s2_patch=slip_model_in(:,13);

SLIP=[s1_patch;s2_patch];
Npatch=length(indx_patch);
dxf=lp/2;
dyf=0;
dzf=0;
theta=(90.0-strkp)*d2r;

[dx,dy]=xy2XY(dxf,dyf,-theta);
dz=dzf;
xxo=xp+dx;
yyo=yp+dy;
zzo=zp+dz;

uE1=0;
uN1=0;
uZ1=0;

uE2=0;
uN2=0;
uZ2=0;

HF=1;

zout=NaN(size(zin));
indx_good=~isnan(zin);
xinsar=XIN(indx_good);
yinsar=YIN(indx_good);
ve=look_e(indx_good);
vn=look_n(indx_good);
vz=look_z(indx_good);


for k=1:Npatch
   xpt=xinsar-xxo(k);
   ypt=yinsar-yyo(k);
   U1=s1_patch(k);
   U2=s2_patch(k);
   delta=dip0(k)*d2r;
   d=-zzo(k);
   len=lp(k);
   W=wp(k);
   fault_type1=1;
   fault_type2=2;
   strike=strkp(k)*d2r;
   tp=zeros(size(xpt));
   [ue1,un1,uz1]=calc_okada(HF,U1,xpt,ypt,nu,delta,d,len,W,fault_type1,strike,tp);
   [ue2,un2,uz2]=calc_okada(HF,U2,xpt,ypt,nu,delta,d,len,W,fault_type2,strike,tp);
   uE1=uE1+ue1;
   uN1=uN1+un1;
   uZ1=uZ1+uz1;

   uE2=uE2+ue2;
   uN2=uN2+un2;
   uZ2=uZ2+uz2;
end

ue=uE1+uE2;
un=uN1+uN2;
uz=uZ1+uZ2;

zout_good=ue.*ve+un.*vn+uz.*vz;
zout(indx_good)=zout_good;

end