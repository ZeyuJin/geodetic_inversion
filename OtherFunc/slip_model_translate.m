function slip_model_out=slip_model_translate(slip_model_in,strk_x,dxf,dyf,dzf);
%
% translate the given slip model by dxf, dyf and dzf
% along strike, dip and vertical directions, respectively
%
% Usage: slip_model_out=slip_model_translate(slip_model_in,strk_x,dxf,dyf,dzf);
%
% slip_model_in:  original (input) fault model
% strk_x: mean strike of the fault
% dxf:  shift along the strike
% dyf:  shift along the dip
% dzf:  shift along vertical %(positive for down shift)
%
% slip_model_out: translated fault model
%
% by Kang Wang on 07/25/2015
% Last updated by Kang Wang on 07/25/2015

format long
d2r=pi/180;

slip_model_out=slip_model_in;
theta_x=(90.0-strk_x)*d2r;

xp_in=slip_model_in(:,4);
yp_in=slip_model_in(:,5);
zp_in=slip_model_in(:,6);

[dxe,dyn]=xy2XY(dxf,dyf,-theta_x);

xp_out=xp_in+dxe;
yp_out=yp_in+dyn;
zp_out=zp_in-dzf;

slip_model_out(:,4)=xp_out;
slip_model_out(:,5)=yp_out;
slip_model_out(:,6)=zp_out;
