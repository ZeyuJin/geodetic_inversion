 function comb_fault = slip2relax(data, filename, slip_mask, fault_index)
% 
% Usage: out=slip2relax(in);
%
% read in a slip model and write it to the format of relax 
%
% note that all the dimensions in relax out will be in meters
% 
% newx and newy are new geographic locations (in meters)
%
% by Kang Wang in June, 1016
% 
% slip_mask means the amplitude that coseismic slip is masked
% usually it is arount 0.5 ~ 1m

d2r=pi/180;

% read original data
indx_fault = data(:,1);
lp_o=data(:,7);
wp_o=data(:,8);
slip1_o=data(:,12);
slip2_o=data(:,13);
area_before = lp_o .* wp_o;
slip_o=sqrt(slip1_o.^2+slip2_o.^2)/1.0e5;

%% mask out slip patches with slip < slip_mask
% just mask out patches on the two parallel NE branches
fault_mask = [];
for i = 1:length(fault_index)
    tmp = data(indx_fault == fault_index(i),:);
    fault_mask = [fault_mask;tmp];
end
slip1_mask = fault_mask(:,12);
slip2_mask = fault_mask(:,13);
slip_all = sqrt(slip1_mask.^2 + slip2_mask.^2);
mask = slip_all >= slip_mask;  % find the bug that masks most coseismic slip
fault_new = fault_mask(mask,:);

% % why substract newx, newy here?
% fault_new(:,4) = fault_new(:,4) - newx;
% fault_new(:,5) = fault_new(:,5) - newy;

% we do not need to mask patches on the rest of fault segments
fault_origin = [];
all_fault = unique(indx_fault);
unmask_index = setdiff(all_fault,fault_index);
for j = 1:length(unmask_index)
    tmp = data(indx_fault == unmask_index(j),:);
    fault_origin = [fault_origin;tmp];
end

% combined fault segments, ignore the number of fault patches index
comb_fault = [fault_new;fault_origin];

%% convert combined fault parameters to RELAX output (km, year, MPa)
xp=comb_fault(:,4);
yp=comb_fault(:,5);
zp=comb_fault(:,6);
lp=comb_fault(:,7);
wp=comb_fault(:,8);
strike_p=comb_fault(:,9);
dip=comb_fault(:,10);
slip1=comb_fault(:,12);
slip2=comb_fault(:,13);

xout=yp/1000;  %north
yout=xp/1000;  %east
zout=-zp/1000; %positive downward

lp_out=lp/1000; 
wp_out=wp/1000;
strike_out=strike_p;
dip_out=dip;
rake_out=atan2(slip2,slip1)/d2r;

% scale the slip amplitude using geodetic moment
slip_out=sqrt(slip1.^2+slip2.^2)/1.0e5;  %convert from cm to km
area_after = lp .* wp;
ratio = sum(slip_o.*area_before) ./ sum(slip_out.*area_after);
slip_out = ratio .* slip_out;

% output for test
disp(ratio);
disp([length(slip_o), length(slip_out)]);

indx_out = [1:length(xout)]';
data_out=[indx_out,slip_out,xout,yout,zout,lp_out,wp_out,strike_out,dip_out,rake_out];

f_out=fopen(filename,'w');
fprintf(f_out,'%5d %12.5e %12.5f %12.5f %12.5f %12.5f %12.5f %8.2f %8.2f %8.2f\n',data_out');
fclose(f_out);

end