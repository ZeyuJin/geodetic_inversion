function [V1,V2,V3,V4]=get_slip_okada_xyz(slip_model_in);
% get the coordinates of patches of a slip model
%
% Usage: [V1,V2,V3,V4]=get_slip_okada_xyz(slip_data);
%
% contents of V1,V2,V3,V4;
%  V1=[x1,y1,z1];
%  V2=[x2,y2,z2];
%  V3=[x3,y3,z3];
%  V4=[x4,y4,z4];
%
% by Kang Wang in Feb. 2016

format long
d2r=pi/180;
ifault=slip_model_in(:,1);
ipatch=slip_model_in(:,2);
idepth=slip_model_in(:,3);
xp=slip_model_in(:,4);
yp=slip_model_in(:,5);
zp=slip_model_in(:,6);
lp=slip_model_in(:,7);
wp=slip_model_in(:,8);
strkp=slip_model_in(:,9);
dip0=slip_model_in(:,10);

Npatch=length(ipatch);
xx1=zeros(Npatch,1);
yy1=zeros(Npatch,1);
zz1=zeros(Npatch,1);

xx2=zeros(Npatch,1);
yy2=zeros(Npatch,1);
zz2=zeros(Npatch,1);

xx3=zeros(Npatch,1);
yy3=zeros(Npatch,1);
zz3=zeros(Npatch,1);

xx4=zeros(Npatch,1);
yy4=zeros(Npatch,1);
zz4=zeros(Npatch,1);

for i=1:Npatch;
  p_theta=(90.0-strkp(i))*d2r;
  p_dip=dip0(i)*d2r;  

  [x1f,y1f]=xy2XY(xp(i),yp(i),p_theta);
    %y1f=y1f+zr(i)*cos(p_dip(i));
   z1f=zp(i);

    x2f=x1f+lp(i);
    y2f=y1f;
    z2f=z1f;

    x3f=x2f;
    y3f=y2f-wp(i)*cos(p_dip);
    z3f=z2f-wp(i)*sin(p_dip);

    x4f=x1f;
    y4f=y3f;
    z4f=z3f;

    [x1,y1]=xy2XY(x1f,y1f,-p_theta);
    [x2,y2]=xy2XY(x2f,y2f,-p_theta);
    [x3,y3]=xy2XY(x3f,y3f,-p_theta);
    [x4,y4]=xy2XY(x4f,y4f,-p_theta);
  
    xx1(i)=x1;
    yy1(i)=y1;
    zz1(i)=z1f;    

    xx2(i)=x2;
    yy2(i)=y2;
    zz2(i)=z2f;

    xx3(i)=x3;
    yy3(i)=y3;
    zz3(i)=z3f;

    xx4(i)=x4;
    yy4(i)=y4;
    zz4(i)=z4f;

end

V1=[xx1,yy1,zz1];
V2=[xx2,yy2,zz2];
V3=[xx3,yy3,zz3];
V4=[xx4,yy4,zz4];
