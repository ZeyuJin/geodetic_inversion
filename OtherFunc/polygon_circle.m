function [xout,yout]=polygon_circle(xin,yin,R);
theta=[0:0.05:2*pi]';


x=R*cos(theta);
y=R*sin(theta);

xout=xin+x;
yout=yin+y;