function cov_mtx_gps=cov_gps(Ngps,sig_gps);
% compute the covariance matrix for the GPS data
%
% by Kang Wang on 12/03/2014
% Last Updated by Kang Wang on 12/03/2014
% Last Updated by Kang Wang on 12/04/2014
%
% Ngps  ------ number of GPS observations
% sig_gps  --- uncertainty of the GPS observation
% W_gps   ---- relative weight of the GPS data

cov_mtx_gps=zeros(Ngps,Ngps);
w_data=1./sig_gps;
sig_sum=sum(w_data);

for i=1:Ngps;
   % A(2) of Fialko 2004
  cov_mtx_gps(i,i)=(1/sig_gps(i))/sig_sum;
end
