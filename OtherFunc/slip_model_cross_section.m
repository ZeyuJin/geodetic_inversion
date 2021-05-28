function [y1,y2,z1,z2,m_out]=slip_model_cross_section(slip_model_in,strk_x);
% get the coordinates (y and z) and total moment of the patches 
% along a profile normal to strike of strk_x;
%
%    Usage: 
%       [y1,y2,z1,z2,m_out]=slip_model_cross_section(slip_model_in,strk_x); 
%
% by Kang Wang on 07/30/2015

format long
d2r=pi/180;

indx_fault=slip_model_in(:,1);
indx_patch=slip_model_in(:,2);
indx_depth=slip_model_in(:,3);

xp=slip_model_in(:,4);
yp=slip_model_in(:,5);
zp=slip_model_in(:,6);

lp=slip_model_in(:,7);
wp=slip_model_in(:,8);
strkp=slip_model_in(:,9);
dip0=slip_model_in(:,10);
tp=slip_model_in(:,11);
slip1=slip_model_in(:,12);
slip2=slip_model_in(:,13);

depth_layer=unique(indx_depth);
Ndepth=length(depth_layer);

y1=zeros(Ndepth,1);
y2=zeros(Ndepth,1);

z1=zeros(Ndepth,1);
z2=zeros(Ndepth,1);

m_out=zeros(Ndepth,1);

theta_x=(90.0-strk_x)*d2r;
[xf,yf]=xy2XY(xp,yp,theta_x);   %rotate the coordiantes to the fault coordinates
zf=zp;

nu=30e9; %shear modulus
for i=1:Ndepth;
   indx_this_layer=find(indx_depth==i);
   xf_this_layer=xf(indx_this_layer);
   yf_this_layer=yf(indx_this_layer);
   zf_this_layer=zf(indx_this_layer);
   lp_this_layer=lp(indx_this_layer);
   wp_this_layer=wp(indx_this_layer);
   dip_this_layer=dip0(indx_this_layer);
   strkp_this_layer=strkp(indx_this_layer);
   slip1_this_layer=slip1(indx_this_layer)/100;  % cm to m
   slip2_this_layer=slip2(indx_this_layer)/100;
   dyf=mean(wp_this_layer)*cos(mean(dip_this_layer)*d2r);
   dzf=mean(wp_this_layer)*sin(mean(dip_this_layer)*d2r);
   y1(i)=mean(yf_this_layer);
   y2(i)=y1(i)-dyf;
   z1(i)=mean(zf_this_layer);
   z2(i)=z1(i)-dzf;
   
   a_this_layer=lp_this_layer.*wp_this_layer; %patch areas
   slip_this_layer=sqrt(slip1_this_layer.^2+slip2_this_layer.^2);
   m_out(i)=sum(nu*slip_this_layer.*a_this_layer);
end



