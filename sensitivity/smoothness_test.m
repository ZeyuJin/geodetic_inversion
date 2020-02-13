close all
clc
clear

dip_change_id = [1 2 3 4 5 7];
H_offset = [-0.6 -0.2 0.5 0.5 0.6 -0.4] * 1e3;
dip_angle = acosd(H_offset ./ 10e3);

slip_model_vs = load_fault_one_plane('fault_geometry/model_A_plane','dip_change_id',dip_change_id,'dip',dip_angle);
% ID = 3;  N_layer_cut = 3;
% [~,nL] = smoo1_each_plane(slip_model_vs);
% patch_this_ID = nL(ID,:);
% total_patch_before = sum(sum(nL(1:ID-1,:),2));
% patch_3_layer_this_ID = 1:sum(patch_this_ID(1:N_layer_cut));
% patch_disappear = total_patch_before + patch_3_layer_this_ID;
% patch_rest_layer_this_ID = max(patch_disappear)+1:sum(sum(nL(1:ID,:),2));
% % recompute the patch number and layer
% slip_model_vs(patch_rest_layer_this_ID,2) = slip_model_vs(patch_rest_layer_this_ID,2) - length(patch_disappear);
% slip_model_vs(patch_rest_layer_this_ID,3) = slip_model_vs(patch_rest_layer_this_ID,3) - N_layer_cut;
% slip_model_vs(patch_disappear,:) = [];

curr_id = max(slip_model_vs(:,1));
slip_model_ds = load_fault_dip_shallow('fault_geometry/model_A_splay','fault_id',curr_id,'dip_change_id',[1:7],'dip_depth',[3.2e3,3.2e3,4.2e3,4.2e3,4.2e3,4.2e3,3e3]);

fault_file = 'fault_geometry/modelA_faults';
segment_file = 'fault_geometry/model_A_seg_connect';     intersect_file = 'fault_geometry/model_A_dip_connect';
iter_step = 3;

SF = [1e-2 2e-2 3e-2 4e-2 5e-2 6e-2 7e-2 8e-2 9e-2 1e-1 2e-1 3e-1 4e-1 5e-1];   % smoothness test
RMS = zeros(size(SF));
Rough = zeros(size(SF));
for jj = 1:length(SF)
    lambda = SF(jj);
    [slip_model,rms,model_roughness] = make_fault_from_insar(slip_model_vs,slip_model_ds,iter_step,'smoothness',lambda, ...
                  'segment_smooth_file',segment_file,'intersect_smooth_file',intersect_file,'shallow_dip_id',[9:13]); 
    RMS(jj) = rms;
    Rough(jj) = model_roughness;
end
L_data = [RMS',Rough'];
save('data_sensitivity/L_curve.mat','L_data');

figure; hold on
plot(L_data(:,1),L_data(:,2),'*-','linewidth',10);
scatter(L_data(:,1),L_data(:,2),500,'r','filled','o');
xlabel('\chi^2');
ylabel('Roughness');
% text(L_data(10,1)+0.002,L_data(10,2)+5,' \lambda = 1','fontsize',30,'fontweight','bold');
% annotation('textarrow',[0.2 0.2],[0.3 0.3]);
set(gca,'fontsize',30,'fontweight','bold');
set(gcf,'PaperPositionMode','auto');