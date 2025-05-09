%% Post Processing of Data: Deramp from residuals
%% All the data will be saved as 'res_deramp_low'
%% Didn't deramp for the ALOS rng since it is toooo noisy.
%% Before that, sub-sample the full-res data ('low') so that the size won't explode (sub_full_res_data)
% 1. Compute Full-Resolution Residuals
% 2. Deramp the data based on residuals

clear
addpath(genpath('/Volumes/T7/Research/MyanmarEQ'));

%% Do the inversion with Maynmar Earthquake Data
%Xiaoyu Zou, 05/05/2025
clear

%% remember to set the input displacement unit into cm. 
%% Preliminary Settings
iint=0;
iter_step2=0;
lonc=96.05;
latc=20.75;
ref_lon=96.05;
model_type = 'okada';
fault_file='fault_file.txt';
segfile='seg_connect.txt';
smoothness=5e-1;
dip_change_id=1;
dip_angle=75;
len_top=1.2e3;%original:1.2e3;
l_ratio=1.2;
w_ratio=1.2;
width=30e3;
depth_start=0;%Everytime you change the depth_start, for the layered model, you have to re-calculate the green's function separately based on the okada slip model.
con=[-1 0 0];

%% Build Fault Geometry
slip_model_vs = load_fault_one_plane(fault_file,'dip_change_id',dip_change_id,'dip',dip_angle, ...
                                     'lonc',lonc,'latc',latc,'ref_lon',ref_lon,'depth_start',depth_start,'len_top',len_top);
slip_model_ds=[];

%% Perform the Inversion
[slip_model,RMS_misfit,model_roughness,insar_model1,insar_model2,insar_model3,insar_model4,rms] = make_fault_from_insar(slip_model_vs,slip_model_ds,iter_step2, ...
                     'segment_smooth_file',segfile,'intersect_smooth_file',[],'fault',fault_file, ...
                     'lonc',lonc,'latc',latc,'ref_lon',ref_lon,'model_type',model_type,'smoothness',smoothness,'con',con);

% compute the residuals
[xo,yo]=ll2xy(lonc,latc,ref_lon);
u=[slip_model(:,12);slip_model(:,13)];

file_path={'SEN/A70/azo_ll_low.grd';'SEN/A143/azo_ll_low.grd';'SEN/D33/azo_ll_low.grd';'SEN/D106/azo_ll_low.grd'};
file_flag=[3,3,3,3]; %1: los; 2: rng; 3: azo
code=95;
for i=4:4
    file=file_path{i};
    [savepath,~,~]=fileparts(file);
    [xx,yy,z]=grdread2(file);
    if file_flag(i)==1
        [~,~,lke]=grdread2(strcat(savepath,"/look_e_low.grd"));
        [~,~,lku]=grdread2(strcat(savepath,"/look_u_low.grd"));
        [~,~,lkn]=grdread2(strcat(savepath,"/look_n_low.grd"));      
    elseif file_flag(i)==2
        [~,~,lke]=grdread2(strcat(savepath,"/look_e_rng_low.grd"));
        [~,~,lku]=grdread2(strcat(savepath,"/look_u_rng_low.grd"));
        [~,~,lkn]=grdread2(strcat(savepath,"/look_n_rng_low.grd"));
    elseif file_flag(i)==3
        [~,~,lke]=grdread2(strcat(savepath,"/look_e_azo_low.grd"));
        [~,~,lku]=grdread2(strcat(savepath,"/look_u_azo_low.grd"));
        [~,~,lkn]=grdread2(strcat(savepath,"/look_n_azo_low.grd"));
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
    if file_flag(i)==3
        G=calc_green_AZO_okada(slip_model,insar_data);
    else
        G=calc_green_insar_okada(slip_model,insar_data);
    end
    insar_residual=Z-G*u;
    insar_residual=reshape(insar_residual,length(yy),length(xx));
    [ramp,cffs]=deramp_xyz(insar_residual,x,y,code);
    if file_flag(i)==1
        grdwrite2(xx,yy,z-ramp,strcat(savepath,'/res_deramped_low.grd'));
    elseif file_flag(i)==2
        grdwrite2(xx,yy,z-ramp,strcat(savepath,'/res_deramped_rng_low.grd'));
    elseif file_flag(i)==3
        grdwrite2(xx,yy,z-ramp,strcat(savepath,'/res_deramped_azo_low.grd'));
    end
end

