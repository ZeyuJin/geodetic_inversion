function w_gps=calc_weight_gps_error(sig_gps);
% calculate the weighting matrix of the GPS data
% Usage: w_gps=calc_weight_gps_error(sig_gps);
%
% by Kang Wang on 08/20/2015
% Last Updated by Kang Wang on 08/20/2015
% Last Updated by Kang Wang on 08/29/2015

format long
Ngps=length(sig_gps);
w_gps=zeros(Ngps,Ngps);
w_data=1./sig_gps;
sig_sum=sum(w_data);

for i=1:Ngps;
    %A(2) of Fialko 2004
    w_gps(i,i)=(1/sig_gps(i))/sig_sum;
end
