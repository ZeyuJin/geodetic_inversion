function [M0,Mw]=slip2moment(slip_model,G);
% Usage: [M0,Mw]=slip2moment(slip_model,G);
%

data_slip=slip_model;
%data_slip=load(slip_model);
indx_patch=data_slip(:,2);
xp=data_slip(:,4);
yp=data_slip(:,5);
zp=data_slip(:,6);
lp=data_slip(:,7);
wp=data_slip(:,8);
strkp=data_slip(:,9);
dip0=data_slip(:,10);
s1_patch=data_slip(:,12)/100;
s2_patch=data_slip(:,13)/100;

A=lp.*wp;
%G=30e9;
D=sqrt(s1_patch.^2+s2_patch.^2);

M0=sum(G*A.*D);
Mw=(2/3)*log10(M0)-6.07;
