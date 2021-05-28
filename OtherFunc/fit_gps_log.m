function [out]=fit_gps_log(tin,uin,sig_in,tau,tout);

ndata=length(tin);
A=zeros(ndata,2);
W=zeros(ndata,2);
B=uin;
% ndata
% size(tin)

for i=1:ndata;
   A(i,1)=1;
   A(i,2)=log(1+tin(i)/tau);
   W(i,i)=1/sig_in(i);
end

C=W*A;
D=W*B;

x=pinv(C)*D;

u_unit=x(1);
amp=x(2);

out=u_unit+amp*log(1+tout/tau);
