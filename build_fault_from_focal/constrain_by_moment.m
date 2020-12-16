function [area,M0] = constrain_by_moment(slip_model,model_type,magnitude,patch_indx)
    fault_id = slip_model(:,1);
    Np = length(fault_id);
    area = zeros(1,2*Np);
    lpatch = slip_model(patch_indx,7);
    wpatch = slip_model(patch_indx,8);    
    Apatch = lpatch' .* wpatch' / 1e2 / 1e6;
    area(patch_indx) = -Apatch;   % slip are negative values    
    M0 = 10^(9.1+1.5*magnitude) / 1e10 / 1e6;
    
    depth = [0 3 6 9 12 15 20 29];
    shear = [2.32 2.71 3.06 3.21 3.35 3.52 3.76 4.25];
    
    if strcmp(model_type,'okada')
        M0 = M0/3;    % the whole shear modulus is 30GPa
    else
        dz_top = slip_model(patch_indx,6) ./ 1000;
        dip = slip_model(patch_indx,10);
        patch_center = dz_top*(-1) + wpatch.*sind(dip)/2/1000;
        for kk = 1:length(patch_center)
            for mm = 2:length(depth)
                if depth(mm) > patch_center(kk)
                    this_depth = depth(mm);
                    top_depth = depth(mm-1);
                    bottom_depth = depth(mm);
                    top_shear = shear(mm-1);
                    bottom_shear = shear(mm);
                    mu_this_depth = (top_shear*(bottom_depth-this_depth)+bottom_shear*(this_depth-top_depth)) ...
                                    /(bottom_depth-top_depth);
                    area(kk) = area(kk) * mu_this_depth;   % scaled by the shear modulus
%                     area(kk+Np) = area(kk+Np) * mu_this_depth;
                    break
                end
            end
        end
    end
end