clear
%getedgrn('edgrnTK')
addpath(genpath('/Volumes/T7/Research/PamirProject'))
% load('slip_model_okada.mat','slip_model');
% 
% load('/Volumes/T7/Research/PamirProject/Inversion/test_data/asc/los_samp0.mat','sampled_insar_data');
% G_raw=calc_green_LOS_edcmp(slip_model,sampled_insar_data,'edgrnPamir');
% save('/Volumes/T7/Research/PamirProject/Inversion/test_data/asc/layered_green.mat','G_raw');
% 
% load('/Volumes/T7/Research/PamirProject/Inversion/test_data/des/los_samp0.mat','sampled_insar_data');
% G_raw=calc_green_LOS_edcmp(slip_model,sampled_insar_data,'edgrnPamir');
% save('/Volumes/T7/Research/PamirProject/Inversion/test_data/des/layered_green.mat','G_raw');


load('slip_model_okada.mat','slip_model');

load('/Volumes/T7/Research/PamirProject/Inversion/test_data/asc/los_samp0.mat','sampled_insar_data');
G_raw=calc_green_LOS_edcmp(slip_model,sampled_insar_data,'edgrnPamir');
save('/Volumes/T7/Research/PamirProject/Inversion/test_data/asc/layered_green.mat','G_raw');

load('/Volumes/T7/Research/PamirProject/Inversion/test_data/des/los_samp0.mat','sampled_insar_data');
G_raw=calc_green_LOS_edcmp(slip_model,sampled_insar_data,'edgrnPamir');
save('/Volumes/T7/Research/PamirProject/Inversion/test_data/des/layered_green.mat','G_raw');

