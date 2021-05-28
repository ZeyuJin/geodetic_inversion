function [SSD_h,depth_patch] = model2slip_profile(slip_model, depth_bins)

    % read model parameters
    zstart = slip_model(:,6);
    lp = slip_model(:,7);
    wp = slip_model(:,8);
    dip = slip_model(:,10);
    uh_strike = slip_model(:,12) ./ 100;  % convert to meters
    uh_dip = slip_model(:,13) ./ 100;  % convert to meters
    uh_total = sqrt(uh_strike.^2 + uh_dip.^2);

    zcenter = zstart - wp/2.*sind(dip);
    zcenter = zcenter * (-1) ./ 1000;
    
    
    % depth_bins
    SSD_h = zeros(length(depth_bins)-1,1);
    depth_patch = zeros(length(depth_bins)-1,1);
    
    for ii = 1:length(depth_bins)-1
        upper_depth = depth_bins(ii);
        bottom_depth = depth_bins(ii+1);
        indx = find(zcenter >= upper_depth & zcenter < bottom_depth);
    
        if isempty(indx)
            SSD_h(ii) = 0;
            depth_patch(ii) = 0;
        else  
            lp_this_depth = lp(indx);
            uh_this_depth = uh_total(indx);
    
            CSP_h = abs(uh_this_depth);   % .* lp_this_depth;
            SSD_h(ii) = sum(CSP_h) / length(CSP_h);
    
            this_depth = zcenter(indx);
            depth_patch(ii) = mean(this_depth);
        
        end
    end
    
    % ignore the data gap in certain depth
    indx_non_zero = find(SSD_h);
    SSD_h = SSD_h(indx_non_zero);
    depth_patch = depth_patch(indx_non_zero);
    
end