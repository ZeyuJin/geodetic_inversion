function compute_moment(slip_model,model_type)
    % compute the scalar seismic moment
    strike_u = slip_model(:,12) ./ 100;     % in meters
    strike_d = slip_model(:,13) ./ 100;
    dz_top = slip_model(:,6) ./ 1000;
    dip = slip_model(:,10);
    D = sqrt(strike_u.^2 + strike_d.^2);
    lpatch = slip_model(:,7);
    wpatch = slip_model(:,8);
    Apatch = lpatch .* wpatch;   
    patch_center = dz_top*(-1) + wpatch.*sind(dip)/2/1000;
    
    depth = [0 3 6 9 12 15 20 29];
    shear = [2.32 2.71 3.06 3.21 3.35 3.52 3.76 4.25] * 10e9; 
    
    if strcmp(model_type,'okada')
        mu = 33e9;
        M0 = sum(mu .* D .* Apatch);
        Mw = 2/3*(log10(M0) - 9.1);
        disp(['The moment magnitude is Mw = ',num2str(Mw)]);
        fprintf('\n');
    else
        total = 0;
        for kk = 1:length(patch_center)
            for mm = 2:length(depth)
                if depth(mm) > patch_center(kk)
                    this_depth = depth(mm);
                    top_depth = depth(mm-1);
                    bottom_depth = depth(mm);
                    top_shear = shear(mm-1);
                    bottom_shear = shear(mm);
                    mu_this_depth = (top_shear*(bottom_depth-this_depth)+bottom_shear*(this_depth-top_depth))/(bottom_depth-top_depth);
                    tmp = D(kk) * Apatch(kk) * mu_this_depth;  
                    break
                end
            end
        total = total + tmp;
        end
        Mw = 2/3*(log10(total) - 9.1);
        disp(['The moment magnitude is Mw = ',num2str(Mw)]);
        fprintf('\n');
    end   
end