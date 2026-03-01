%% Save the Forward Models for Sampling
%% Before that, sub-sample the full-res data ('low') so that the size won't explode (sub_full_res_data), except for ALOS RNG
%% The output is 'model.grd'


clear
addpath(genpath('/Volumes/T7/Research/Tingri_project'));

%% 1. Compute the Inversion
iint=0;
iter_step2=0;
lonc=87.5;
latc=28.7;
ref_lon=87.5;

load slip_model.mat
% compute the residuals
[xo,yo]=ll2xy(lonc,latc,ref_lon);
%u=[slip_model(:,2);slip_model(:,3)];
% file_path={'SEN/des121/data/los_ll_masked_deramped.grd','SEN/des48/data/los_ll_masked_deramped.grd','SEN/asc12/data/los_ll_deramped.grd','ALOS/data/los_ll_masked_deramped.grd',...
%     'SEN/des121/data/rng/rng_ll_masked.grd','SEN/asc12/data/rng/rng_ll_masked.grd'};
file_path={'SEN/des121/data/res_deramped_low.grd','SEN/des48/data/res_deramped_low.grd','SEN/asc12/data/res_deramped_low.grd','ALOS/data/res_deramped_low.grd',...
    'SEN/des121/data/rng/res_deramped_low.grd','SEN/asc12/data/rng/res_deramped_low.grd','ALOS/data/rng/rng_ll_low.grd'};
file_flag=[1,1,1,1,2,2,2]; %1: los; 2: rng; 

for i=3:7
    file=file_path{i};
    [savepath,~,~]=fileparts(file);
    [xx,yy,z]=grdread2(file);
    if file_flag(i)==1
        [~,~,lke]=grdread2(strcat(savepath,"/look_e_low.grd"));
        [~,~,lku]=grdread2(strcat(savepath,"/look_u_low.grd"));
        [~,~,lkn]=grdread2(strcat(savepath,"/look_n_low.grd"));      
    else
        [~,~,lke]=grdread2(strcat(savepath,"/look_rng_e_low.grd"));
        [~,~,lku]=grdread2(strcat(savepath,"/look_rng_u_low.grd"));
        [~,~,lkn]=grdread2(strcat(savepath,"/look_rng_n_low.grd"));
    end
    [X,Y]=meshgrid(xx,yy);
    x=X;
    y=Y;
    X=X(:);
    Y=Y(:);
    [X,Y]=ll2xy(X,Y,ref_lon);
    X=X-xo;
    Y=Y-yo;
    Z=z(:);
    lke=lke(:);
    lku=lku(:);
    lkn=lkn(:);
    insar_data=[X,Y,Z,lke,lkn,lku];
    G=calc_green_insar_okada(slip_model,insar_data);
    insar_model=G*u;
    insar_model=reshape(insar_model,length(yy),length(xx));
    grdwrite2(xx,yy,insar_model,strcat(savepath,'/model.grd'));
end
