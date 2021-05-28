function hh = scatter_patches(varargin)
%scatter_patches Display scatter plot using patch objects 
%
%   This function supports all Name/Value properties supported by patch 
%   objects for scatter markers.  This may be especially useful for 
%   setting the 'FaceAlpha' property, enabling alpha-blended scatter plots.
%
%   scatter_patches(X,Y,S,C) displays colored circles at the locations 
%   specified by the vectors X and Y (which must be the same size).  
%
%   S determines the area of each marker (in points^2). S can be a
%   vector the same length a X and Y or a scalar. If S is a scalar, 
%   we draw all the markers the same size. If S is empty, the
%   default size is used.
%   
%   C determines the colors of the markers. When C is a vector the
%   same length as X and Y, the values in C are linearly mapped
%   to the colors in the current colormap. When C is a 
%   length(X)-by-3 matrix, it directly specifies the colors of the  
%   markers as RGB values. C can also be a color string. See ColorSpec.
%
%   scatter_patches(X,Y) draws the markers in the default size and color.
%
%   scatter_patches(X,Y,S) draws the markers at the specified sizes (S)
%   with a single color. This type of graph is also known as
%   a bubble plot.
%   scatter_patches(...,M) uses marker M. Available markers:
%      o    Circle (default)
%      d    Diamond
%      s    Square
%      ^    Triangle (pointing up)
%      v    Triangle (pointing down)
%      <    Triangle (pointing left)
%      >    Triangle (pointing right)
%
%   scatter_patches(X,Y,...,'Name','Value','Name','Value'...)
%     passes the Name/Value argument to the patch objects.  For example:
%       scatter_patches(randn(100,1),randn(100,1),'r','FaceAlpha',0.3);
%
%   scatter_patches(AX,...) plots into AX instead of GCA.
%
%   H = scatter_patches(...) returns handles to the patch objects created.
%
%   Example:
%      N=100;
%      hh1=scatter_patches(randn(N,1),1*randn(N,1),36, 'r','FaceAlpha',0.4,'EdgeColor','none');
%      hold on;
%      hh2=scatter_patches(randn(N,1),2*randn(N,1),36, N*rand(N,1),'s','FaceAlpha',0.2,'EdgeColor','none');
%      hold on;
%      hh3=scatter_patches(randn(N,1),4 + 1*randn(N,1),100*rand(N,1), 1:N,'<','FaceAlpha',0.2,'EdgeColor','r');
%      legend([hh1(1),hh2(1),hh3(1)], {'red circles','multicolor squares','multicolor triangles'});
%
% Artemy Kolchinsky, Indiana University, 2014


[cax, args] = axescheck(varargin{:});
if isempty(cax) || ishghandle(cax,'axes')
    cax = newplot(cax);
else
    cax = ancestor(cax,'Axes');
end

[args,pvpairs] = parseparams(args);

nargs = length(args);
if nargs < 2
    error('Not enough arguments supplied (need at least x,y)');
end

xs = args{1};
ys = args{2};

if any([numel(xs) numel(ys) length(xs) length(ys)] ~= length(xs))
    error('X and Y must be single-dimensional vectors of same length');
end

if nargs > 2
    sizes = args{3};
    if numel(sizes) ~= 1 && numel(sizes) ~= numel(xs)
        error('Marker size parameter must be single entry or vector with same number of elements as points');
    end
    if any(sizes < 0)
        error('Sizes must be positive');
    end
else
    sizes = 36;
end

sizes = sqrt(sizes);

colors = [];
if nargs > 3
    colors = args{4};
else
    if ~isempty(pvpairs) && length(pvpairs{1})==1
        colors = pvpairs{1};
        pvpairs = pvpairs(2:end);
    end
end

if nargs > 4
    error('Too many arguments');
end


if isempty(colors)
    colors = 'b';
end
if numel(colors) ~= 1 && numel(colors) ~= length(xs) && numel(colors) ~= 3 && numel(colors) ~= 3*length(xs)
    error('Colors  parameter must be single entry or vector with same number of elements as points, or RGB triplet or Nx3 matrix of RGB ');
end
if numel(colors) == length(xs) && size(colors,1) == 1,
    colors = colors';
end


markerType = 'o';
if ~isempty(pvpairs) && length(pvpairs{1})==1
    markerType = pvpairs{1};
    pvpairs = pvpairs(2:end);
end

if mod(numel(pvpairs),2) ~= 0
    error('Need an even number of Name,Value parameters');
end

patch_rotation = 0;
if markerType == 'o'
    patch_num_points = 30;
elseif markerType == 's'
    patch_num_points = 4;
    patch_rotation = pi/4;
elseif markerType == 'd'
    patch_num_points = 4;
elseif markerType == '<'
    patch_num_points = 3;
    patch_rotation = 3*pi/2;
elseif markerType == '>'
    patch_num_points = 3;
    patch_rotation = pi/2;
elseif markerType == '^'
    patch_num_points = 3;
elseif markerType == 'v'
    patch_num_points = 3;
    patch_rotation = pi;
else
    error(['Unknown markerType ' markerType]);
end

patchSpec = (0:(2*pi/patch_num_points):2*pi) + patch_rotation;

baseArgs = {};
perPointArgs = {};
for j=1:2:numel(pvpairs)
    if size(pvpairs{j+1},1)==1
        baseArgs{end+1} = pvpairs{j};
        baseArgs{end+1} = pvpairs{j+1};
    else
        if size(pvpairs{j+1},1) ~= length(xs)
            error(sprintf('Argument %s must have single entry or the same number of rows as number of points', pvpairs{j}));
        end
        perPointArgs{end+1} = pvpairs{j};
        perPointArgs{end+1} = pvpairs{j+1};
    end
end

hh=[];
cUnits = get(cax,'Units');
set(cax,'Units','points');
pos = get(cax,'Position');
xlims = get(cax,'XLim');
ylims = get(cax,'YLim');
ptsPerXUnit = pos(3)/diff(xlims);
ptsPerYUnit = pos(4)/diff(ylims);
for i=1:numel(xs)
    cPatchArgs = baseArgs;
    for j=1:2:numel(perPointArgs)
        cPatchArgs{end+1} = perPointArgs{j};
        cPatchArgs{end+1} = perPointArgs{j+1}(i,:);
    end
    if length(sizes) > 1
        cSize = sizes(i);
    else
        cSize = sizes;
    end
    if size(colors,1) > 1
        cColor = colors(i,:);
    else
        cColor = colors;
    end
    hh(end+1) = patch( cSize * sin(patchSpec) / ptsPerXUnit + xs(i), cSize * cos(patchSpec) / ptsPerYUnit + ys(i), cColor, cPatchArgs{:});
end

set(cax,'Units',cUnits);

callBackFcn = @(varargin) updatePatches(cax,hh,xs,ys,sizes,patchSpec);
hl = addlistener(handle(cax),{'XLim','YLim','Position','DataAspectRatio'},'PostSet',callBackFcn);
set(get(cax,'Parent'),'ResizeFcn',callBackFcn);

set(cax,'UserData', {get(cax,'UserData') hl});


function updatePatches(cax, patchObjects, xs, ys, sizes, patchSpec)
    if ~ishandle(patchObjects(1))
        % already deleted
        return
    end
    cUnits = get(cax,'Units');
    set(cax,'Units','points');
    pos = get(cax,'Position');
    xlims = get(cax,'XLim');
    ylims = get(cax,'YLim');
    ptsPerXUnit = pos(3)/diff(xlims);
    ptsPerYUnit = pos(4)/diff(ylims);
    for i=1:numel(xs)
        if length(sizes) > 1
            cSize = sizes(i);
        else
            cSize = sizes;
        end
        set(patchObjects(i), 'XData', cSize * sin(patchSpec) / ptsPerXUnit + xs(i));
        set(patchObjects(i), 'YData', cSize * cos(patchSpec) / ptsPerYUnit + ys(i));
    end

    set(cax,'Units',cUnits);
