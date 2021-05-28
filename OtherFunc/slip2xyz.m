function [UX,UY,UZ]=slip2xyz(xpt,ypt,zpt,slip_model);
%function to calculate the surface displacements at given points using slip
%model with Okada's solutions;
%****************INPUT**************
%xpt  ---- vector of  x-coordinates of observations
%ypt  ---- vector of y-coordinates of observations
%zpt --- vector of z-coordindates of observations

%slip_model ---- file name of the slip model
% Last Updated by Kang Wang on 01/21/2015

d2r=pi/180;

%show_slip_model(slip_model);  %show slip model
data_slip=load(slip_model);
indx_patch=data_slip(:,2);
xp=data_slip(:,4);
yp=data_slip(:,5);
zp=data_slip(:,6);

lp=data_slip(:,7);
wp=data_slip(:,8);
strkp=data_slip(:,9);
dip0=data_slip(:,10);
s1_patch=data_slip(:,12);   %slip along strike
s2_patch=data_slip(:,13);   % thrust component
spatch=sqrt(s1_patch.^2+s2_patch.^2);  %slip magnitude
rake_patch=atan2(s2_patch,s1_patch)/d2r;

Npatch=length(indx_patch);
dxf=lp/2;
dyf=-cos(dip0*d2r).*wp/2;
dzf=-sin(dip0*d2r).*wp/2;

theta=(90.0-strkp)*d2r;
[dx,dy]=xy2XY(dxf,dyf,-theta);

xo=xp+dx;
yo=yp+dy;     %coordinates of the center of the patches
zo=zp+dzf;

Npt=length(xpt);
UX=zeros(Npt,1);
UY=zeros(Npt,1);
UZ=zeros(Npt,1);

for i=1:Npatch;
    XE=xpt-xo(i);
    YN=ypt-yo(i);
    Zpatch=-zo(i);
    [uE,uN,uZ]=okada85(XE,YN,Zpatch,strkp(i),dip0(i),lp(i),wp(i),rake_patch(i),spatch(i),0,0.25);
    UX=UX+uE;
    UY=UY+uN;
    UZ=UZ+uZ;
end
