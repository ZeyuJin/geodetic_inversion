function build_fault_M5(lon_eq,lat_eq,depth,strike,dip,Mw,varargin)
% build only one segment based on the focal mechanism of USGS
% lon_eq,lat_eq,depth is the origin of the fault (not to be the same as the hypocenter)
% strike,dip,Mw are from the first-motion solution (from USGS)
% apply a gaussian taper to the slip distribution

    % hypocenter from the focal mechanism
    % lon_eq = -117.575;
    % lat_eq = 35.760;
    % depth = 7e3;
    % strike = 218.68;
    % dip = 68;
    
    % default values
    lonf = -117.5;   % reference point
    latf = 35.5;
    [xc,yc] = utm2ll(lon_eq,lat_eq,0,1);
    len = 7e3;
    W = 7e3;
    Npatch = 15;
    sigma_l = 5e3;
    sigma_w = 5e3;  
    cmax = 30;
    model_name = 'M5_slip_gauss';
    model_type = 'okada';
    slip_type = 'right';
    axis_range = [-10 5 17 30 -15 0];
    seis_matfile = [];
    fault_file = [];
    
    % read varargin values and assembly
    if ~isempty(varargin)
        for CC = 1:floor(length(varargin)/2)
            try
                switch lower(varargin{CC*2-1})
                    case 'lonf'
                        lonf = varargin{CC*2};
                    case 'latf'
                        latf = varargin{CC*2};
                    case 'fault_length'
                        len = varargin{CC*2};
                    case 'fault_width'
                        W = varargin{CC*2};
                    case 'npatch'
                        Npatch = varargin{CC*2};  % Npatch*Npatch rectangulars
                        % better to define the size Lpatch*Wpatch
                        % in the future
                    case 'sig_len'
                        sigma_l = varargin{CC*2};
                    case 'sig_width'
                        sigma_w = varargin{CC*2};
                    case 'color_range'
                        cmax = varargin{CC*2};
                    case 'save_name'
                        model_name = varargin{CC*2};
                    case 'model_type'
                        model_type = varargin{CC*2};
                    case 'slip_type'
                        slip_type = varargin{CC*2};
                    case 'axis_range'
                        axis_range = varargin{CC*2};
                    case 'seismic' 
                        seis_matfile = varargin{CC*2};
                    case 'fault_file'
                        fault_file = varargin{CC*2};
                end
            catch
                error('Unrecognized Keyword');
            end
        end
    end
    
    [xo,yo] = utm2ll(lonf,latf,0,1);    % reference point
    fxc = xc - xo;
    fyc = yc - yo;

    % build one segment centered on the M5 hypocenter
    xf1 = fxc + len/2*sind(strike-180);
    yf1 = fyc + len/2*cosd(strike-180);
%     xf2 = fxc - len/2*sind(strike-180);
%     yf2 = fyc - len/2*cosd(strike-180);

    % define the fault segment geometry 7x7 km locally (default)
    A = 1;
    len_inv = len/Npatch;
    W_inv = W/Npatch;
    xp = len_inv/2*[1:2:Npatch*2-1] - len/2;  % origin at the center patch
    yp = W_inv/2*[1:2:Npatch*2-1] - W/2;
    [xpm,ypm] = meshgrid(xp,yp);

    % apply the 2-D Gaussian distribution centered on the middle patch
    % at the centers of all patches
    xo = 0; yo = 0;    % middle of the patch
    G = A * exp(-(xpm-xo).^2/2/sigma_l^2 - (ypm-yo).^2/2/sigma_w^2);

    % multiply the Gaussian by a constant F such that the sum of the geodetic moment 
    % over all patches equals moment that corresponds to Mw5.4.
    patch_area = W*len/Npatch^2;
    M0 = 10^(1.5*Mw+9.1);
%     F = M0 / mu / sum(patch_area*G(:));
%     slip_dist = F*G;

    % % test for the plot
    % figure;
    % pcolor(xpm,ypm,F*G);
    % shading flat
    % colormap jet
    % colorbar

    % save in the model.mat file
    xstart = xf1 - W/2*cosd(dip)*cosd(360-strike);     % fix the center bug
    ystart = yf1 - W/2*cosd(dip)*sind(360-strike);
    dz_start = -depth+W/2*sind(dip);
    strike_input = strike;
    dip_input = dip;
    N_layer = Npatch;
    lp_top = len_inv;
    bias_lp = 1;     % uniform patches
    bias_wp = 1;
    fault_id = 1;
    fault_M5 = make_fault_segments(fault_id,xstart,ystart,dz_start,strike_input,dip_input, ...
                              len,W,N_layer,lp_top,bias_lp,bias_wp);
    
    if strcmp(model_type,'okada')
        mu = 30e9;
        F = M0 / mu / sum(patch_area.*G(:));    
    else
        M5_indx = 1:size(fault_M5,1);
        [layer_area,M0] = constrain_by_moment(fault_M5,model_type,Mw,M5_indx);
        if strcmp(slip_type,'left')
            layer_area = layer_area * -1;    % left lateral slip is positive
        end
        Np = length(layer_area) / 2;
        patch_area = layer_area(1:Np);    % assume strike-slip only
        % each column indicates the slip of each layer, because the fault
        % models starts from each layer
        G_tmp = G';  
        F = M0 / sum(patch_area'.*G_tmp(:)) / 100;
    end
    
    slip_dist = F*G;
    fault_M5(:,2) = [1:size(fault_M5,1)]';
    fault_M5(:,12) = 100 * reshape(slip_dist',1,[]);  % assume strike-slip only
    slip_model = fault_M5;
    
    % test for plotting the segment
    show_slip_model(fault_M5,'misfit_range',cmax,'axis_range',axis_range,'seismic',seis_matfile,'fault_file',fault_file);    
    
    read_plot_fault_segment('../verified_rupture_trace.txt');
    read_plot_fault_segment('../not_verified_rupture_trace.txt');
    
    title('M5 event');
    save([model_name,'.mat'],'slip_model');

end