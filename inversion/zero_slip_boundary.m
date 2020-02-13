function [W,d] = zero_slip_boundary(slip_model,segment_ID,top_layer_no,ratio)
% add zero slip constraint at the top surface
    Np = size(slip_model,1);
    all_fault_id = slip_model(:,1);
    all_patch_id = slip_model(:,2);
    all_layer_id = slip_model(:,3);
    
%     N_patch_before = length(find(all_fault_id<segment_ID));
    patch_this_segment = all_patch_id(all_fault_id==segment_ID);
    layer_this_segment = all_layer_id(all_fault_id==segment_ID);
    patch_top_this_segment = patch_this_segment(layer_this_segment <= top_layer_no);
    
    strike_indx = patch_top_this_segment';
%     disp(strike_indx);
    dip_indx = strike_indx + Np;
%     disp(dip_indx);
    zero_slip_indx = [strike_indx,dip_indx];
    
    W = zeros(2*Np);
    for ii = 1:length(zero_slip_indx)
        W(zero_slip_indx(ii),zero_slip_indx(ii)) = ratio;
    end
    d = zeros(2*Np,1);
end