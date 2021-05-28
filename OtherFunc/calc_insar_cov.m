function covd = calc_insar_cov(xin,yin,sigma,L);
% calculate the covariacne matrix for a given set of observations
% assuming that the c(r) = sigma* exp(-r/L);
%
% Usage: covd = calculate_insar_cov(xin,yin,sigma,L);
%
% xin ---- x coordinates of the observations
% yin ---- y coordinates of the observations
%
% by Kang Wang in Feb. 2018

npt=length(xin);
covd=zeros(npt,npt);

for k=1:npt;
 xnow = xin(k);
 ynow = yin(k);
 
 dx=xnow -xin;
 dy=ynow -yin;
 dr=sqrt(dx.^2+dy.^2);
 covd(k,:)=sigma*exp(-dr/L);
end