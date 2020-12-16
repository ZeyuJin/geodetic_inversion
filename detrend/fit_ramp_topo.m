function [pout] = fit_ramp_topo(p,x,y,h);
% fit the input phase as:
%
%  phs = a*x + b*y +c*topo+d;
%  p=[a,b,c,d];
%
% Usage :  p=fit_ramp_topo(phs,x,y,h);
%  

indx_good=~isnan(p);
phi=p(indx_good);
Nin=length(phi);
xin=x(indx_good);
yin=y(indx_good);
hgt=h(indx_good);

G=zeros(Nin,4);

for k=1:Nin;
 G(k,1)=xin(k);
 G(k,2)=yin(k);
 G(k,3)=hgt(k);
 G(k,4)=1;
end

mm=pinv(G)*phi;
a1=mm(1);
b1=mm(2);
c1=mm(3);
d1=mm(4);

pout=[a1,b1,c1,d1];
