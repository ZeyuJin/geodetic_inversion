function [pout] = fit_ramp_topo(p,x,y,h)
% fit the input phase as:
%
%  phs = a*x + b*y +c*topo+d;
%  p=[a,b,c,d];
%
% Usage :  p=fit_ramp_topo(phs,x,y,h);
%  

indx_good=~isnan(p);
phi=p(indx_good);
xin=x(indx_good);
yin=y(indx_good);
hgt=h(indx_good);

indx2=~isinf(phi);
phi=phi(indx2);
xin=xin(indx2);
yin=yin(indx2);
hgt=hgt(indx2);

% in case some NaNs in DEM
indx3=~isnan(hgt);
phi=phi(indx3);
xin=xin(indx3);
yin=yin(indx3);
hgt=hgt(indx3);

Nin = length(phi);
G=zeros(Nin,4);

% for k=1:Nin
 G(:,1)=xin(:);
 G(:,2)=yin(:);
 G(:,3)=hgt(:);
 G(:,4)=1;
% end

mm=pinv(G)*phi;
% mm = G\phi;   % speed up instead of SVD
a1=mm(1);
b1=mm(2);
c1=mm(3);
d1=mm(4);

pout=[a1,b1,c1,d1];

end