function slip_model_out=slip_model_rotate_dip(slip_model_in,strk_x,dyf_rotate,d_angle);
% rotate a slip model along the axis going through yf_rotate by d_angle
%
% Usage: slip_model_out=slip_model_rotate_dip(slip_model_in,strk_x,yf_rotate,d_angle;
% slip_model_in: horizontal slip model
% yf_rotate: horizontal distance from the origin along  dip direciton 
% d_angle: downward rotation [deg]
%
% by Kang Wang on 07/27/2015
% Last Updated by Kang Wang on 07/27/2015
% Last Updated by Kang Wang on 10/25/2016

format long
d2r=pi/180;

theta_x=(90.0-strk_x)*d2r;
xp_in=slip_model_in(:,4);
yp_in=slip_model_in(:,5);
zp_in=slip_model_in(:,6);
[xf_in,yf_in]=xy2XY(xp_in,yp_in,theta_x);
yf_in=-yf_in+dyf_rotate;
%yf_in=-yf_in;

%zf_in=zp_in+0e3;
zf_in=0e3;

dtheta=d_angle*d2r;
[yf_out,zf_out]=xy2XY(yf_in,zf_in,dtheta);
xf_out=xf_in;
yf_out=-yf_out+dyf_rotate;

[xp_out,yp_out]=xy2XY(xf_out,yf_out,-theta_x);
slip_model_out=slip_model_in;
slip_model_out(:,4)=xp_out;
slip_model_out(:,5)=yp_out;
slip_model_out(:,6)=zf_out+zp_in;
slip_model_out(:,10)=d_angle;
