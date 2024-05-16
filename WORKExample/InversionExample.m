%% Geodetic Inversion Test
%11/15/2022, Xiaoyu Zou
%For Turkey earthquake inversion, try to start from Step 4
clear
clc
%% Step 0: set path and read in parameters & files
addpath(genpath('/Users/x3zou/Desktop/Turkey_test/geodetic_inversion-master'))
setenv('PATH',[getenv('PATH'),':/opt/homebrew/bin']);  % add the path of GMT
configfile='configfile.txt';
fid = fopen(configfile);
tmp_txt = fgetl(fid);
files=zeros(1,1);
files=string(files);
while tmp_txt ~=-1
    tmp_txt = fgetl(fid);
    skip=find(tmp_txt=='#');
    if ~isempty(skip)
        continue
    end
    files=[files tmp_txt];
end
files(1)=[];
files(end)=[];
fclose(fid)
grdin=files(1);
grdout=files(2);
data_list=files(3);  % LOS data
los_list=data_list;
fault_file=files(4);
segment_smooth_file = files(5);
segment_file = files(6);
intersect_smooth_file = [];
intersect_file=[];
slip_model_ds=[];%set empty for most cases, because it's used to construct the geometry of "Y shape" or "flower structure" formed by shallow splay faults.



configpara='configpara.txt';
fid = fopen(configpara);
tmp_txt = fgetl(fid);
para=zeros(1,1);
para=string(para);
while tmp_txt ~=-1
    tmp_txt = fgetl(fid);
    skip=find(tmp_txt=='#');
    if ~isempty(skip)
        continue
    end
    para=[para tmp_txt];
end
para(1)=[];
para(end)=[];
para=str2double(para);
fclose(fid)
lonf=para(1);
latf=para(2);%coordinate of pixel point
ref_lon=para(3);
threshold=para(4);
lonc=para(5);
latc=para(6);%coordinate of reference point (0,0)
Nmin=para(7);
Nmax=para(8);
dip_change_id=para(9):para(10);%the array of fault ids that have dip angles not equal to 90 degrees
dip_angle=[para(11) para(12) para(13) para(14) para(15)];%the array of fault ids that have dip angles not equal to 90 degrees
iter_step=para(16);
iter_step2=para(17);



%% Step 1: data cleaning using clean_insar-data. Remove some near-field
%%unwrapping erros manually first.
clean_insar_data

%% Step 2: detrend the phase and remove the phase ambiguity
grdin1='/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/Research/TurkeyEQ/WORKExample/DES5/LOS2/unwrap_clean_sample.grd';
grdout1='/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/Research/TurkeyEQ/WORKExample/DES5/LOS2/los_clean_detrend.grd';
remove_ref_from_grid(grdin1,grdout1,lonf,latf,ref_lon,threshold);
grdin2='/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/Research/TurkeyEQ/WORKExample/ASC100/LOS2/unwrap_clean_sample.grd';
grdout2='/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/Research/TurkeyEQ/WORKExample/ASC100/LOS2/los_clean_detrend.grd';
remove_ref_from_grid(grdin2,grdout2,lonf,latf,ref_lon,threshold);

%% Step 3: apply quad-tree sampling to all detrended data 
make_insar_data(data_list,Nmin,Nmax,'method','quadtree','fault',fault_file,'area',[71.8 73.9 37.7 39.1],'ref_lon',71);


%% Step 4: Build the fault geometry
slip_model_vs = load_fault_one_plane(fault_file,'dip_change_id',dip_change_id,'dip',dip_angle,'lonc',72,'latc',38.5,'ref_lon',71,'len_top',1.2e3);


%% Step 5: inversion using first downsampled data
[slip_model,~,~,insar_model1,insar_model2] = make_fault_from_insar3(slip_model_vs,slip_model_ds,iter_step,'shallow_dip_id',[],'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'fault',fault_file,'lonc',72,'latc',38.5,'ref_lon',71,'model_type','okada');
iint=iter_step;
plot_insar_model_resampled(['/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/Research/TurkeyEQ/WORKExample/ASC100/LOS2/los_samp',num2str(iint),'.mat'],insar_model1,'iter_step',iint,'fault',fault_file,'model_type','okada','misfit_range',30,'defo_max',120,'ref_lon',ref_lon,'lonc',-117.5,'latc',35.5);
plot_insar_model_resampled(['/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/Research/TurkeyEQ/WORKExample/DES5/LOS2/los_samp',num2str(iint),'.mat'],insar_model2,'iter_step',iint,'fault',fault_file,'model_type','okada','misfit_range',30,'defo_max',120,'ref_lon',ref_lon,'lonc',-117.5,'latc',35.5);
show_slip_model(slip_model,'misfit_range',400,'ref_lon',ref_lon,'lonc',-117.5,'latc',35.5,'axis_range',[60 110 -45 25 -25 0]);


%% Step 6: iterative sampling data using the model predictions (Wang and Fialko, GRL 2015) 
resamp_insar_data(slip_model,data_list,Nmin, Nmax, iter_step2, 'fault', fault_file, 'dec',2, 'lonc',72, 'latc',38.5, 'ref_lon',71);

%% Step7: inversion using resampled data
iint=iter_step2;
[slip_model,rms,model_roughness,insar_model1,insar_model2] = make_fault_from_insar3(slip_model_vs,slip_model_ds,iter_step2, ...
                     'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'fault',fault_file, ...
                     'lonc',72,'latc',38.5,'ref_lon',71);
plot_insar_model_resampled(['/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/Research/TurkeyEQ/WORKExample/ASC100/LOS2/los_samp',num2str(iint),'.mat'],insar_model1,'iter_step',iint,'fault',fault_file,'model_type','okada','misfit_range',30,'defo_max',120,'ref_lon',ref_lon,'lonc',-117.5,'latc',35.5);
plot_insar_model_resampled(['/Users/xiaoyuzou/Library/CloudStorage/OneDrive-UCSanDiego/Research/TurkeyEQ/WORKExample/DES5/LOS2/los_samp',num2str(iint),'.mat'],insar_model2,'iter_step',iint,'fault',fault_file,'model_type','okada','misfit_range',30,'defo_max',120,'ref_lon',ref_lon,'lonc',-117.5,'latc',35.5);
show_slip_model(slip_model,'misfit_range',400,'ref_lon',ref_lon,'lonc',-117.5,'latc',35.5,'axis_range',[60 110 -45 25 -25 0]);

