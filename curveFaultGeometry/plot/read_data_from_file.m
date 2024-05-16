function [ux, uy, uz] = read_data_from_file(filename)
    addpath('/Users/zej011/work/Kang_tutorial/codes_utilities/matlab/igppsar');
    
    ref_lon = 99;
    lonc = 99;
    latc = 34;
    [xo, yo] = ll2xy(lonc, latc, ref_lon);
    
    data = load(filename);
    lon = data(:,1);
    lat = data(:,2);
    [ux, uy] = ll2xy(lon, lat, ref_lon);
    ux = (ux - xo) ./ 1e3;
    uy = (uy - yo) ./ 1e3;
    
    ncols = size(data, 2);
    if ncols >= 3
        uz = data(:,3) * (-1);
    else
        uz = zeros(size(ux));  % for surface traces
    end
    
end