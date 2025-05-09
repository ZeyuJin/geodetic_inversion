function layered_pscmp(filename,greens_dir,out_dir,snapshot_file,nrec,O_lat,O_lon,O_depth,len,W,strike,dip,slps,slpd,varargin)
% solve the bug that dip slip component should be reversed 
% when converting the hanging wall to the footwall 
% Because EDCMP did not allow dip angles larger than 90 degrees
% that is: strike_after = strike_before - 180
%          dip_after = 180 - dip_before
%          dip_slip_after = -1 * dip_slip_before

obs_type = 0;      % default is irregular observation points
obs_arr = [];      % observation array boundaries, in format of [xr1,yr1,xr2,yr2];
if ~isempty(varargin)
    for CC = 1:floor(length(varargin)/2)
        try
            switch lower(varargin{CC*2-1})
                case 'obs_type'
                    obs_type = varargin{CC*2};
                case 'obs_arr'
                    obs_arr = varargin{CC*2};
            end
        catch
            error('Unrecognized Keyword');
        end
    end
end

fedc = fopen(filename,'wt');
fprintf(fedc,['#===============================================================================\n' ...
'# OBSERVATION ARRAY\n' ...
'# =================\n' ...
'# 1. selection for irregular observation positions (= 0) or a 1D observation\n' ...
'#    profile (= 1) or a rectangular 2D observation array (= 2): iposrec\n' ...
'#\n' ...
'#    IF (iposrec = 0 for irregular observation positions) THEN\n' ...
'#\n' ...
'# 2. number of positions: nrec\n' ...
'#\n' ...
'# 3. coordinates of the observations: (lat(i),lon(i)), i=1,nrec\n' ...
'#\n' ...
'#    ELSE IF (iposrec = 1 for regular 1D observation array) THEN\n' ...
'#\n' ...
'# 2. number of position samples of the profile: nrec\n' ...
'#\n' ...
'# 3. the start and end positions: (lat1,lon1), (lat2,lon2)\n' ...
'#\n' ...
'#    ELSE IF (iposrec = 2 for rectanglular 2D observation array) THEN\n' ...
'#\n' ...
'# 2. number of x samples, start and end values: nxrec, xrec1, xrec2\n' ...
'#\n' ...
'# 3. number of y samples, start and end values: nyrec, yrec1, yrec2\n' ...
'#\n' ...
'#    sequence of the positions in output data: lat(1),lon(1); ...; lat(nx),lon(1);\n' ...
'#    lat(1),lon(2); ...; lat(nx),lon(2); ...; lat(1),lon(ny); ...; lat(nx),lon(ny).\n' ...
'#\n' ...
'#    Note that the total number of observation positions (nrec or nxrec*nyrec)\n' ...
'#    should be <= NRECMAX (see pecglob.h)!\n' ...
'#===============================================================================\n']);
fprintf(fedc,' %d \n',obs_type);
if obs_type == 0 || obs_type == 1
    fprintf(fedc,' %d \n',nrec);
    if obs_type == 1
        fprintf(fedc,'(%.4f  %.4f), (%.4f  %.4f)\n',obs_arr(3:4),obs_arr(1:2));
    end
elseif obs_type == 2
    if length(nrec) ~= 2, error('Something wrong with receiver points!'); end
    fprintf(fedc,' %d  %.4f  %.4f\n',nrec(2),obs_arr(2),obs_arr(4));
    fprintf(fedc,' %d  %.4f  %.4f\n',nrec(1),obs_arr(1),obs_arr(3));
else 
    error('There is something wrong with observation type');
end


fprintf(fedc,['#===============================================================================\n' ...
'# OUTPUTS\n' ...
'#\n' ...
'# 1. select (1/0) output for los displacement (only for snapshots, see below),\n' ...
'#    x, y, and z-cosines to the INSAR orbit: insar, xlos, ylos, zlos \n' ...
'#\n' ...
'#    if this option is selected, the snapshots will include additional data: \n' ...
'#    LOS_Dsp = los displacement to the given satellite orbit.\n' ...
'#\n' ...
'# 2. select (1/0) output for Coulomb stress changes (only for snapshots, see \n' ...
'#    below): icmb, friction, Skempton ratio, strike, dip, and rake angles [deg] \n' ...
'#    describing the uniform regional master fault mechanism, the uniform regional \n' ...
'#    principal stresses: sigma1, sigma2 and sigma3 [Pa] in arbitrary order (the \n' ...
'#    orietation of the pre-stress field will be derived by assuming that the\n' ...
'#    master fault is optimally oriented according to Coulomb failure criterion) \n' ...
'#\n' ...
'#    if this option is selected (icmb = 1), the snapshots will include additional \n' ...
'#    data: \n' ...
'#    CMB_Fix, Sig_Fix = Coulomb and normal stress changes on master fault; \n' ...
'#    CMB_Op1/2, Sig_Op1/2 = Coulomb and normal stress changes on the two optimally \n' ...
'#                       oriented faults; \n' ...
'#    Str_Op1/2, Dip_Op1/2, Slp_Op1/2 = strike, dip and rake angles of the two \n' ...
'#                       optimally oriented faults. \n' ...
'# \n' ...
'#    Note: the 1. optimally orieted fault is the one closest to the master fault. \n' ...
'# \n' ... 
'# 3. output directory in char format: outdir \n' ...
'# \n' ...
'# 4. select outputs for displacement components (1/0 = yes/no): itout(i), i=1,3 \n' ...
'# \n' ...
'# 5. the file names in char format for the x, y, and z components: \n' ...
'#    toutfile(i), i=1,3 \n' ...
'# \n' ...
'# 6. select outputs for stress components (1/0 = yes/no): itout(i), i=4,9 \n' ...
'# \n' ...
'# 7. the file names in char format for the xx, yy, zz, xy, yz, and zx components: \n' ...
'#    toutfile(i), i=4,9 \n' ...
'# \n' ...
'# 8. select outputs for vertical NS and EW tilt components, block rotation, geoid \n' ...
'#    and gravity changes (1/0 = yes/no): itout(i), i=10,14 \n' ...
'# \n' ...
'# 9. the file names in char format for the NS tilt (positive if borehole top \n' ...
'#    tilts to north), EW tilt (positive if borehole top tilts to east), block \n' ...
'#    rotation (clockwise positive), geoid and gravity changes: toutfile(i), i=10,14 \n' ...
'# \n' ...
'#    Note that all above outputs are time series with the time window as same \n' ...
'#    as used for the Green''s functions \n' ...
'# \n' ...
'#10. number of scenario outputs ("snapshots": spatial distribution of all above \n' ...
'#    observables at given time points; <= NSCENMAX (see pscglob.h): nsc \n' ...
'# \n' ...
'#11. the time [day], and file name (in char format) for the 1. snapshot; \n' ...
'#12. the time [day], and file name (in char format) for the 2. snapshot; \n' ...
'#13. ... \n' ...
'# \n' ...
'#    Note that all file or directory names should not be longer than 80 \n' ...
'#    characters. Directories must be ended by / (unix) or \\ (dos)! \n' ...
'#===============================================================================\n']);
fprintf(fedc,' 0 \n');
fprintf(fedc,' 0 \n');
% system(['mkdir -p ',out_dir]);
fprintf(fedc,['''',out_dir,'/''\n']);
fprintf(fedc,' 0             0             0 \n');
fprintf(fedc,' ''ux.dat''    ''uy.dat''    ''uz.dat'' \n');
fprintf(fedc,' 0             0             0             0              0             0 \n');
fprintf(fedc,' ''sxx.dat''   ''syy.dat''   ''szz.dat''   ''sxy.dat''    ''syz.dat''   ''szx.dat'' \n');
fprintf(fedc,' 0             0             0             0              0\n');
fprintf(fedc,' ''tx.dat''    ''ty.dat''    ''rot.dat''   ''gd.dat''     ''gr.dat'' \n');

% read the snapshots file
fid = fopen(snapshot_file,'r');
CC = textscan(fid,'%.2f %s %s\n');
day_time = CC{1};
snapshot_name = CC{2};
% comment = CC{3};
fclose(fid);

ndays = length(day_time);
fprintf(fedc,' %d \n',ndays);
cell_time = num2cell(day_time);
mycell = {cell_time{:};snapshot_name{:}};
fprintf(fedc,'  %.2f    %s   \n',mycell{:});
fprintf(fedc,['#===============================================================================\n' ...
'#\n' ...
'# GREEN''S FUNCTION DATABASE \n' ...
'# ========================= \n' ...
'# 1. directory where the Green''s functions are stored: grndir \n' ...
'#\n' ...
'# 2. file names (without extensions!) for the 13 Green''s functions: \n' ...
'#    3 displacement komponents (uz, ur, ut): green(i), i=1,3 \n' ...
'#    6 stress components (szz, srr, stt, szr, srt, stz): green(i), i=4,9 \n' ...
'#    radial and tangential components measured by a borehole tiltmeter, \n' ...
'#    rigid rotation around z-axis, geoid and gravity changes (tr, tt, rot, gd, gr): \n' ...
'#    green(i), i=10,14 \n' ...
'# \n' ...
'#    Note that all file or directory names should not be longer than 80 \n' ...
'#    characters. Directories must be ended by / (unix) or \\ (dos)! The \n' ...
'#    extensions of the file names will be automatically considered. They \n' ...
'#    are ".ep", ".ss", ".ds" and ".cl" denoting the explosion (inflation) \n' ...
'#    strike-slip, the dip-slip and the compensated linear vector dipole \n' ...
'#    sources, respectively. \n' ...
'#\n' ...
'#===============================================================================\n']);
fprintf(fedc,['''',greens_dir,'/''\n']);
fprintf(fedc,' ''uz'' ''ur'' ''ut'' \n');
fprintf(fedc,' ''szz'' ''srr'' ''stt'' ''szr'' ''srt'' ''stz'' \n');
fprintf(fedc,' ''tr''  ''tt''  ''rot'' ''gd''  ''gr'' \n');

fprintf(fedc,['#===============================================================================\n' ...
'# RECTANGULAR SUBFAULTS \n' ...
'# ===================== \n' ...
'# 1. number of subfaults (<= NSMAX in pscglob.h), latitude [deg] and east \n' ...
'#    longitude [deg] of the regional reference point as  origin of the Cartesian \n' ...
'#    coordinate system: ns, lat0, lon0 \n' ...
'#\n' ...
'# 2. parameters for the 1. rectangular subfault: geographic coordinates \n' ...
'#    (O_lat, O_lon) [deg] and O_depth [km] of the local reference point on \n' ...
'#    the present fault plane, length (along strike) [km] and width (along down \n' ...
'#    dip) [km], strike [deg], dip [deg], number of equi-size fault \n' ...
'#    patches along the strike (np_st) and along the dip (np_di) (total number of \n' ...
'#    fault patches = np_st x np_di), and the start time of the rupture; the \n' ...
'#    following data lines describe the slip distribution on the present sub- \n' ...
'#    fault: \n' ...
'#\n' ...
'#    pos_s[km]  pos_d[km]  slip_along_strike[m]  slip_along_dip[m]  opening[m] \n' ...
'#\n' ...
'#    where (pos_s,pos_d) defines the position of the center of each patch in \n' ...
'#    the local coordinate system with the origin at the reference point: \n' ...
'#    pos_s = distance along the length (positive in the strike direction) \n' ...
'#    pos_d = distance along the width (positive in the down-dip direction)\n' ...
'#\n' ...
'# 3. ... for the 2. subfault ...\n' ...
'# ...\n' ...
'#                   N \n' ...
'#                  / \n' ...
'#                 /| strike \n' ...
'#                +-------------------------\n' ...
'#                |\\        p .            \\ W \n' ...
'#                :-\\      i .              \\ i \n' ...
'#                |  \\    l .                \\ d \n' ...
'#                :90 \\  S .                  \\ t \n' ...
'#                |-dip\\  .                    \\ h \n' ...
'#                :     \\. | rake               \\ \n' ...
'#                Z      --------------------------\n' ...
'#                              L e n g t h \n' ...
'#\n' ...
'#    Note that a point inflation can be simulated by three point openning \n' ...
'#    faults (each causes a third part of the volume of the point inflation) \n' ...
'#    with orientation orthogonal to each other. the results obtained should \n' ...
'#    be multiplied by a scaling factor 3(1-nu)/(1+nu), where nu is the Poisson \n' ...
'#    ratio at the source. The scaling factor is the ratio of the seismic \n' ...
'#    moment (energy) of an inflation source to that of a tensile source inducing \n' ...
'#    a plate openning with the same volume change. \n' ...
'#===============================================================================\n' ...
'# n_faults (Slip model by Ji Chen, USGS) \n' ...
'#-------------------------------------------------------------------------------\n']);

nfaults = length(strike);
fprintf(fedc,'  %d \n',nfaults);
fprintf(fedc,['#-------------------------------------------------------------------------------\n' ...
'# n   O_lat   O_lon    O_depth length  width strike dip   np_st np_di start_time \n' ...
'# [-] [deg]   [deg]    [km]    [km]     [km] [deg]  [deg] [-]   [-]   [day] \n' ...
'#     pos_s   pos_d    slp_stk slp_dip open \n' ...
'#     [km]    [km]     [m]     [m]     [m] \n' ...
'#-------------------------------------------------------------------------------\n']);
pos_s = len/2;
pos_d = W/2;
for kk = 1:nfaults
    % PSCMP does not allow dip angles larger than 90 degrees
    % also strike angles not larger than 360 degrees
    % change the hanging wall to the corresponding footwall
    if dip(kk) > 90
        this_strk = strike(kk) - 180;
        this_dip = 180 - dip(kk);
        this_slpd = slpd(kk) * -1;  % dip-slip is relative with each other
        % while strike-slip has absolute direction in local coordinates.
    else
        this_strk = strike(kk);
        this_dip = dip(kk);
        this_slpd = slpd(kk);
    end
    fprintf(fedc,' %d  %.4f  %.4f  %.4f  %.4f  %.4f  %.2f  %.2f  1  1  0.00\n', ...
                 kk,O_lat(kk),O_lon(kk),O_depth(kk),len(kk),W(kk),this_strk,this_dip);
    fprintf(fedc,'     %.4f  %.4f  %.5f  %.5f  0.00\n',pos_s(kk),pos_d(kk),slps(kk),this_slpd);
end

fprintf(fedc,'#================================end of input===================================');
status = fclose(fedc);
if status == 0
    disp('The file is successfully written!');
end
end
