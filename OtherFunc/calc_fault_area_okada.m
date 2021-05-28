function S=calc_fault_area_okada(slip_model_in);
%  Calculate the total fault area from the given slip model
% 
%  Usage: S=calc_fault_area_okada(slip_model_in);
%
%  by Kang Wang on 08/25/2015

format long
lp=slip_model_in(:,7);
wp=slip_model_in(:,8);
sp=lp.*wp;
S=sum(sp);  %in the unit of m^2;
