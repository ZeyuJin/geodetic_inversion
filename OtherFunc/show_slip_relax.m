function sho_slip_relax(fault_data)
d2r=pi/180;

ipatch=fault_data(:,1);
%slip=fault_data(:,2);
%xp=fault_data(:,4);
%yp=fault_data(:,3);
%zp=-fault_data(:,5);
%lp=fault_data(:,6);
%wp=fault_data(:,7);
%strkp=fault_data(:,8);
%dip=fault_data(:,9);
%rake=fault_data(:,10);

n=length(ipatch);
hold on
for i=1:n;
  slip=fault_data(i,2);
  xp=fault_data(i,4);
  yp=fault_data(i,3);
  zp=-fault_data(i,5);
  lp=fault_data(i,6);
  wp=fault_data(i,7);
  strike=fault_data(i,8);
  dip=fault_data(i,9);
  rake=fault_data(i,10);
    p_theta=(90.0-strike)*d2r;

   [x1f,y1f]=xy2XY(xp,yp,p_theta);
   z1f=zp;

   x2f=x1f+lp;
   y2f=y1f;
   z2f=z1f;

   x3f=x2f;
   y3f=y2f-wp*cos(dip*d2r);
   z3f=z2f-wp*sin(dip*d2r);

   x4f=x1f;
   y4f=y3f;
   z4f=z3f;

   [x1,y1]=xy2XY(x1f,y1f,-p_theta);
    [x2,y2]=xy2XY(x2f,y2f,-p_theta);
    [x3,y3]=xy2XY(x3f,y3f,-p_theta);
    [x4,y4]=xy2XY(x4f,y4f,-p_theta);

    X=[x1,x2,x3,x4];
    Y=[y1,y2,y3,y4];
    Z=[z1f,z2f,z3f,z4f];

    C=[slip,slip,slip,slip];
    patch(X,Y,Z,C);
end
%axis equal
colorbar
colormap('jet');
hold off
end
