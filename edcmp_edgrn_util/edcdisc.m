function [pxs,pys,pzs,pmoment]=edcdisc(xs_patch,ys_patch,zs_patch,slip_patch,lp_patch, ...
    wp_patch,strike_patch,dip_patch,rake_patch,nz,z1,z2,dr,dz)

    % discretize different point sources along each patch
    % strike and dip angle are degrees as the input
    d2r = pi / 180;
    sm = zeros(3,3);

    st = strike_patch * d2r;
    di = dip_patch * d2r;
    ra = rake_patch * d2r;

    if di > 0
        dwidth = min(dr, dz/sin(di));
    else
        dwidth = dr;
    end
    dlength = dr;

    % coefficients of moment tensor
    sm(1,1) = -sin(di)*cos(ra)*sin(2*st) - sin(2*di)*sin(ra)*sin(st)^2;
    sm(2,2) = sin(di)*cos(ra)*sin(2*st) - sin(2*di)*sin(ra)*cos(st)^2;
    sm(3,3) = -(sm(1,1) + sm(2,2));
    sm(1,2) = sin(di)*cos(ra)*cos(2*st) + 0.5*sin(2*di)*sin(ra)*sin(2*st);
    sm(2,1) = sm(1,2);
    sm(1,3) = -cos(di)*cos(ra)*cos(st) - cos(2*di)*sin(ra)*sin(st);
    sm(3,1) = sm(1,3);
    sm(2,3) = -cos(di)*cos(ra)*sin(st) + cos(2*di)*sin(ra)*cos(st);
    sm(3,2) = sm(2,3);

    nx = max(1, round(lp_patch/dlength));
    ny = max(1, round(wp_patch/dwidth));
    dx = lp_patch / nx;
    dy = wp_patch / ny;

    % if both length and width = 0, then it is a point source
    disarea = slip_patch;
    if dx > 0
        disarea = disarea * dx;
    end
    if dy > 0
        disarea = disarea * dy;
    end

    % initialize empty array
    pmoment = zeros(5,nx*ny);
    % all points are located exactly within each line/rectangles
    % It's better not to locate any points at the edges of mesh
    x = dx * ((1:nx) - 0.5);
    y = dy * ((1:ny) - 0.5);
    [x, y] = meshgrid(x, y);
    x = x(:);
    y = y(:);
    pxs = xs_patch + x*cos(st) - y*cos(di)*sin(st);
    pys = ys_patch + x*sin(st) + y*cos(di)*cos(st);
    pzs = zs_patch + y*sin(di);
    
    ind_shallow = find(pzs < z1-dz, 1);
    ind_deeper = find(pzs > z2+dz, 1);
    
    if ~isempty(ind_shallow)
        disp('Warning: parts of source rectangles shallower than the Green function grids!');
    end
    
    if ~isempty(ind_deeper)
        disp('Warning: parts of source rectangles deeper than the Green function grids!');
    end
    
    pmoment(1,:) = sm(1,2) * disarea;
    pmoment(2,:) = sm(1,3) * disarea;
    pmoment(3,:) = sm(3,3) * disarea;
    pmoment(4,:) = 0.5*(sm(1,1)-sm(2,2))*disarea;
    pmoment(5,:) = sm(2,3) * disarea;

%     % can be accelerated using meshgrid
%     nps = 0;
%     for ix = 1:nx
%         x = dx * (ix - 0.5);
%         for iy = 1:ny
%             y = dy * (iy - 0.5);
%             nps = nps + 1;
%             
%             pxs(nps) = xs_patch + x*cos(st) - y*cos(di)*sin(st);
%             pys(nps) = ys_patch + x*sin(st) + y*cos(di)*cos(st);
%             pzs(nps) = zs_patch + y*sin(di);
% 
%             if pzs(nps) < z1 - dz
%                 disp('Warning: parts of source rectangles shallower than the Green function grids!');
%             end
% 
%             if pzs(nps) > z2 + dz
%                 disp('Warning: parts of source rectangles deeper than the Green function grids!');
%             end
% 
%             pmoment(1,nps) = sm(1,2) * disarea;
%             pmoment(2,nps) = sm(1,3) * disarea;
%             pmoment(3,nps) = sm(3,3) * disarea;
%             pmoment(4,nps) = 0.5*(sm(1,1)-sm(2,2))*disarea;
%             pmoment(5,nps) = sm(2,3) * disarea;
%         end
%     end
%     disp(['Total number of point sources for this patch: ', num2str(nps)]);
end