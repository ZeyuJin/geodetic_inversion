function [xout,yout,zout]=grd_detrend(x_in,y_in,z_in);
[X,Y]=meshgrid(x_in,y_in);
xout=x_in;
yout=y_in;

indx_good=~isnan(z_in);
z_good=z_in(indx_good);
x_good=X(indx_good);
y_good=Y(indx_good);

A=[x_good,y_good,ones(size(x_good))];
b=pinv(A)*z_good;

%z_fit=A*b;
zfit=X*b(1)+Y*b(2)+b(3);
zout=z_in-zfit;
