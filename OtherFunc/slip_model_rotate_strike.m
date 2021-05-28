function slip_model_out=slip_model_rotate_strike(slip_model_in,strk_x,dstrk);
%
% rotate the given slip model counter-closewise by dstrk
%
% Usage: slip_model_out=slip_model_rotate_strike(slip_model_in,strk_x,dstrk);
%
% slip_model_in: orginal (input) fault model
% strk_x: mean strke of the original fault model
% dstrk:  degree of counter-closewise ratation (deg)
%
% slip_model_out: updated fault model
%
%by Kang Wang on 07/26/2015
%Last Updated by Kang Wang on 07/26/2015

format long
d2r=pi/180;

slip_model_out=slip_model_in;
theta_x=(90.0-strk_x)*d2r;

xp_in=slip_model_in(:,4);
yp_in=slip_model_in(:,5);
zp_in=slip_model_in(:,6);
strk_in=slip_model_in(:,9);


dtheta=dstrk*d2r;

theta_new=dtheta+theta_x;

[xf_in,yf_in]=xy2XY(xp_in,yp_in,theta_x);

[xp_out,yp_out]=xy2XY(xf_in,yf_in,-theta_new);
zp_out=zp_in;


strk_out=strk_in-dstrk;
slip_model_out(:,4)=xp_out;
slip_model_out(:,5)=yp_out;
slip_model_out(:,6)=zp_out;
slip_model_out(:,9)=strk_out;
