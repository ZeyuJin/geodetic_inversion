function dis=p2p_dis(p1,p2);
% calculate the distance between two points p1(x1,y1,z1), p2(x2,y2,z2);
% by Kang Wang: 03/21/2013
x1=p1(1);
y1=p1(2);
z1=p1(3);

x2=p2(1);
y2=p2(2);
z2=p2(3);

dis=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);

