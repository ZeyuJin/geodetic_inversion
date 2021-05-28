function w_insar=calc_weight_insar_error(rms_insar);
% calculate the weight  matrix of the InSAR from the uncertainity of the data
%  
%  Usage: w_insar=calc_weight_insar_error(rms_insar);
% by Kang Wang on 08/29/2015

format long
Nobs=length(rms_insar);
w_insar=zeros(Nobs,Nobs);
w_data=1./rms_insar;
sig_sum=sum(w_data);

for i=1:Nobs;
   %A(2) of Fialko 2004
   w_insar(i,i)=(1/rms_insar(i))/sig_sum;
end


