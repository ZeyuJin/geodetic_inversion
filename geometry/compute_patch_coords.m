function [x_center,y_center,z_center] = compute_patch_coords(xp,yp,zp,lp,wp,strike,dip)
% the input arguments means the X,Y,Z coordinates of one patch corner
% X,Y start from the strike direction, Z start from the top
% lp,wp mean the length and width of this patch
% the function will compute the central coordinates of this patch
    d2r = pi / 180;
    p_theta = (90 - strike) * d2r;
    p_dip = dip * d2r;
    [x1f,y1f] = xy2XY(xp,yp,p_theta);
    z1f = zp;
    
    x2f = x1f + lp;
    y2f = y1f - wp .* cos(p_dip);
    z2f = z1f - wp .* sin(p_dip);
    
    xco = mean([x1f,x2f]);
    yco = mean([y1f,y2f]);
    zco = mean([z1f,z2f]);
    
    [x_center,y_center] = xy2XY(xco,yco,-p_theta);
    z_center = zco;
end