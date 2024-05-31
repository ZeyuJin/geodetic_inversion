function [Vx, Vy, Vd] = depth_dependent_point(meshGrid, pointGrid, biasL, biasW, W, ratio)
    % used to interpolate depth-dependent points using the meshgrid
    % and the point grid
    % size(meshgrid,1) = size(pointgrid,1) - 1
    % size(meshgrid,2) = size(pointgrid,2) - 1
    % ratio controls the top width of triangle mesh, in order to 
    % make sure at least one point is located between each triangle
    % return interpolated 3D scattered points

    % destruct the grid data
    x1m = meshGrid(:,1);
    d1m = meshGrid(:,2);
    y1m = meshGrid(:,3);
    
    x1p = pointGrid(:,1);
    d1p = pointGrid(:,2);
    y1p = pointGrid(:,3);
 
    topN = numel(unique(x1p));        % top layer of point sources
    d1p_uni = unique(d1p);       
    topW = (d1p_uni(2) - d1p_uni(1)) * ratio;   % top width along dip direction

    % determine how many layers first
    Nest = ceil(log(1 + W/topW*(biasW-1)) / log(biasW));  % estimated number of layer
    wp_factor = [biasW.^(0:Nest-1)];
    W_layer = [0, W ./ sum(wp_factor) .* wp_factor];

    % then determine how many points per layer
    lp_factor = [biasL.^(1:Nest+1)];
    pointLayer = ceil(topN ./ lp_factor);

    
    % compute the position of depth-dependent points
    % even distance along the strike
    % make sure all point sources are located within the convex hull
    Vx = [];  
    Vd = [];
    minX = min(x1m);
    maxX = max(x1m);

    for i = 1:Nest+1
        xin = linspace(minX, maxX, pointLayer(i));

        if i == Nest+1
            tmpD = -W - eps;  % make sure point sources within the grid
        else
            tmpD = -sum(W_layer(1:i));
        end

        din = tmpD .* ones(size(xin));

        Vx = [Vx; xin'];
        Vd = [Vd; din'];
    end

    % interpolation using the dense grid data
    % make sure each vertex is located within the curved surface
    Vy = griddata([x1m; x1p], [d1m; d1p], [y1m; y1p], Vx, Vd, 'cubic');
end
