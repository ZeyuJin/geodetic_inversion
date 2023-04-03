function [xb, yb, zb, xp, yp, zp] = gridFitInterpolate(xfit, yfit, zfit, xin, yin)
    % [xb, yb, zb] are interpolated boundary points
    % [xp, yp, zp] are interpolated center point source locations
    % return values are matrices 
    
    zb = gridfit(xfit, yfit, zfit, xin, yin);
    [xb, yb] = meshgrid(xin, yin);
    
    xp = movmean(xin, 2, 'Endpoints', 'discard');
    yp = movmean(yin, 2, 'Endpoints', 'discard');
    zp = gridfit(xfit, yfit, zfit, xp, yp);
    [xp, yp] = meshgrid(xp, yp);
end