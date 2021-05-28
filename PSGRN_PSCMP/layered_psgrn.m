function layered_psgrn(filename,data_dir,nrec,r1,r2,ndepth,zs1,zs2,time,layer_profile,varargin)				       				       
% configure file for the module PSGRN

samp_rat = 1;   % uniform sampling
wavenumber = 0.025;
grav_eff = 0;
nt = 1;  % only consider the coseismic changes
% read varargin values and assembly
if ~isempty(varargin)
    for CC = 1:floor(length(varargin)/2)
        try
            switch lower(varargin{CC*2-1})
                case 'samp_ratio'
                    samp_rat = varargin{CC*2};
                case 'wavenumber'
                    wavenumber = varargin{CC*2};
                case 'grav_effect'
                    grav_eff = varargin{CC*2};
                case 'time_samples'
                    nt = varargin{CC*2};
            end
        catch
            error('Unrecognized Keyword');
        end
    end
end

fedc = fopen(filename,'wt');
fprintf(fedc,['#------------------------------------------------------------------------------\n' ...
'#       PARAMETERS FOR SOURCE-OBSERVATION CONFIGURATIONS\n' ...
'#       ================================================\n' ...
'# 1. the uniform depth of the observation points [km], switch for oceanic (0)\n' ...
'#    or continental(1) earthquakes;\n' ...
'# 2. number of (horizontal) observation distances (> 1 and <= nrmax defined in\n' ...
'#    psgglob.h), start and end distances [km], ratio (>= 1.0) between max. and\n' ...
'#    min. sampling interval (1.0 for equidistant sampling);\n' ...
'# 3. number of equidistant source depths (>= 1 and <= nzsmax defined in\n' ...
'#    psgglob.h), start and end source depths [km];\n' ...
'#\n' ...
'#    r1,r2 = minimum and maximum horizontal source-observation\n' ...
'#            distances (r2 > r1).\n' ...
'#    zs1,zs2 = minimum and maximum source depths (zs2 >= zs1 > 0).\n' ...
'#\n' ...
'#    Note that the same sampling rates dr_min and dzs will be used later by the\n' ...
'#    program "pscmp08" for discretizing the finite source planes to a 2D grid\n' ...
'#    of point sources.\n' ...
'#------------------------------------------------------------------------------\n']);
fprintf(fedc,'      0.0  1 \n');   % surface observations and continental(1) earthquakes
fprintf(fedc,' %d  %.2f  %.2f  %.2f \n',nrec,r1,r2,samp_rat);  % ratio between max. and min. sampling interval
fprintf(fedc,' %d  %.2f  %.2f \n',ndepth,zs1,zs2);
% for i=1:nrec
%     fprintf(fedc,' ( %e , %e ) \n',yr(i),xr(i));  % x is North in edcmp
% end 

% binary observation points file to be read into PSCMP
% file_rec = [data_dir,'.rec'];

%%%%% write receiver location file (observation points)
% temp = [yr(:), xr(:)]';   % x is North in edcmp
% fid  = fopen([dpath,'/',data_dir,'/',file_rec],'w');
% fwrite(fid, temp, 'real*4');
% fclose(fid);
%%%%% included in each dataset directory

fprintf(fedc,['#------------------------------------------------------------------------------\n' ...
'#       PARAMETERS FOR TIME SAMPLING\n' ...
'#       ============================\n' ...
'# 1. number of time samples (<= ntmax def. in psgglob.h) and time window [days].\n' ...
'#\n' ...
'#    Note that nt (> 0) should be power of 2 (the fft-rule). If nt = 1, the\n' ...
'#    coseismic (t = 0) changes will be computed; If nt = 2, the coseismic\n' ...
'#    (t = 0) and steady-state (t -> infinity) changes will be computed;\n' ...
'#    Otherwise, time series for the given time samples will be computed.\n' ...
'#\n' ...
'#------------------------------------------------------------------------------\n']);
fprintf(fedc,' %d  %.1f \n',nt,time);   % time in days


fprintf(fedc,['#===============================================================================\n' ...
'#       PARAMETERS FOR WAVENUMBER INTEGRATION\n' ...
'#       =====================================\n' ...
'# 1. relative accuracy of the wave-number integration (suggested: 0.1 - 0.01)\n' ...
'# 2. factor (> 0 and < 1) for including influence of earth''s gravity on the\n' ...
'#    deformation field (e.g. 0/1 = without / with 100%% gravity effect).\n' ...
'#------------------------------------------------------------------------------\n']);
fprintf(fedc,' %.3f \n',wavenumber);
fprintf(fedc,' %.3f \n',grav_eff);   % NOT consider gravity effect


fprintf(fedc,['#------------------------------------------------------------------------------\n' ...
'#       PARAMETERS FOR OUTPUT FILES\n' ...
'#       ===========================\n' ...
'#\n' ...
'# 1. output directory\n' ...
'# 2. file names for 3 displacement components (uz, ur, ut)\n' ...
'# 3. file names for 6 stress components (szz, srr, stt, szr, srt, stz)\n' ...
'# 4. file names for radial and tangential tilt components (as measured by a\n' ...
'#    borehole tiltmeter), rigid rotation of horizontal plane, geoid and gravity\n' ...
'#    changes (tr, tt, rot, gd, gr)\n' ...
'#\n' ...
'#    Note that all file or directory names should not be longer than 80\n' ...
'#    characters. Directory and subdirectoy names must be separated and ended\n' ...
'#    by / (unix) or \\ (dos)! All file names should be given without extensions\n' ...
'#    that will be appended automatically by ".ep" for the explosion (inflation)\n' ...
'#    source, ".ss" for the strike-slip source, ".ds" for the dip-slip source,\n' ...
'#    and ".cl" for the compensated linear vector dipole source)\n' ...
'#\n' ...
'#------------------------------------------------------------------------------\n']);
system(['mkdir -p ',data_dir]);
outdir=['''',data_dir,'/''']; 
fprintf(fedc,' %s \n',outdir);
fprintf(fedc,' ''uz'' ''ur'' ''ut'' \n');
fprintf(fedc,' ''szz'' ''srr'' ''stt'' ''szr'' ''srt'' ''stz'' \n');
fprintf(fedc,' ''tr''  ''tt''  ''rot'' ''gd''  ''gr'' \n');


fprintf(fedc,['#------------------------------------------------------------------------------\n' ...
'#       GLOBAL MODEL PARAMETERS\n' ...
'#       =======================\n' ...
'# 1. number of data lines of the layered model (<= lmax as defined in psgglob.h)\n' ...
'#\n' ...
'#    The surface and the upper boundary of the half-space as well as the\n' ...
'#    interfaces at which the viscoelastic parameters are continuous, are all\n' ...
'#    defined by a single data line; All other interfaces, at which the\n' ...
'#    viscoelastic parameters are discontinuous, are all defined by two data\n' ...
'#    lines (upper-side and lower-side values). This input format could also be\n' ...
'#    used for a graphic plot of the layered model. Layers which have different\n' ...
'#    parameter values at top and bottom, will be treated as layers with a\n' ...
'#    constant gradient, and will be discretised to a number of homogeneous\n' ...
'#    sublayers. Errors due to the discretisation are limited within about 5%%\n' ...
'#    (changeable, see psgglob.h).\n' ...
'#\n' ...
'# 2.... parameters of the multilayered model\n' ...
'#\n' ...
'#    Burgers rheology [a Kelvin-Voigt body (mu1, eta1) and a Maxwell body\n' ...
'#    (mu2, eta2) in series connection] for relaxation of shear modulus is\n' ...
'#    implemented. No relaxation of compressional modulus is considered.\n' ...
'#\n' ...
'#    eta1  = transient viscosity (dashpot of the Kelvin-Voigt body; <= 0 means\n' ...
'#            infinity value)\n' ...
'#    eta2  = steady-state viscosity (dashpot of the Maxwell body; <= 0 means\n' ...
'#            infinity value)\n' ...
'#    alpha = ratio between the effective and the unrelaxed shear modulus\n' ...
'#            = mu1/(mu1+mu2) (> 0 and <= 1) (unrelaxed modulus mu2 is\n' ...
'#            derived from S wave velocity and density)\n' ...
'#\n' ...
'#    Special cases:\n' ...
'#        (1) Elastic: eta1 and eta2 <= 0 (i.e. infinity); alpha meaningless\n' ...
'#        (2) Maxwell body: eta1 <= 0 (i.e. eta1 = infinity)\n' ...
'#                          or alpha = 1 (i.e. mu1 = infinity)\n' ...
'#        (3) Standard-Linear-Solid: eta2 <= 0 (i.e. infinity)\n' ...
'#------------------------------------------------------------------------------\n']);
layer_para = load(layer_profile);
n_layer = size(layer_para,1);
ID = layer_para(:,1);  depth = layer_para(:,2);
Vp = layer_para(:,3);  Vs = layer_para(:,4);  rho = layer_para(:,5);
eta1 = layer_para(:,6);  eta2 = layer_para(:,7);  alpha = layer_para(:,8);
fprintf(fedc,' %d \n',n_layer);
fprintf(fedc,'#------------------------------------------------------------------------------\n');
fprintf(fedc,'# no  depth[km]  vp[km/s]  vs[km/s]  rho[kg/m^3] eta1[Pa*s] eta2[Pa*s] alpha\n');
fprintf(fedc,'#------------------------------------------------------------------------------\n');
formatSpec = ' %d  %.3f  %.3f  %.3f  %.1f  %8.2e  %8.2e  %.3f\n';
write_info = [ID,depth,Vp,Vs,rho,eta1,eta2,alpha]';
fprintf(fedc,formatSpec,write_info);
fprintf(fedc,'#=======================end of input===========================================');

status = fclose(fedc);
if status == 0
    disp('The file is successfully written!');
end

end