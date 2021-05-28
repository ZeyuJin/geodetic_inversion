function slip2stress_edcmp(slip_model,save_dir,green_dir,depth)
% this script will compute 6 stress components using EDCMP
% based on the full resolution grid data 
   
   % load the fault geometry and compute the slip and rake direction
   xs = slip_model(:,4);
   ys = slip_model(:,5);
   zs = slip_model(:,6);
   lp = slip_model(:,7);
   wp = slip_model(:,8);
   strike = slip_model(:,9);
   dip = slip_model(:,10);                                                                                                
   U_strike = slip_model(:,12);                                                
   U_dip = slip_model(:,13);
   U = sqrt(U_strike.^2 + U_dip.^2) ./ 100;   % slip in m using EDCMP
   rake = atan2d(U_dip,U_strike);
   
   % use EDCMP to compute the surface 3D displacements
   % format same as DIS3D
   cmdstr = ['mkdir -p ',save_dir,'/outdata'];  % mkdir for the edcmp output
   system(cmdstr);      
   count = '_M5';
   out_type = 'strs';
   obs_type = 2;    % grid observation
%    obs_type = 1;    % line observation
%    obs_arr = [-35 49.5 24.5 -10]*1e3;     % start from the upper left point 
%    obs_arr = [40 30 120 -50]*1e3;
%    obs_arr = [40 150 190 -50]*1e3;
   obs_arr = [0 40 100 -60]*1e3;  % for 2016 M5.4 event
   inv = 1000;       % interval is 500m
   xrec = abs(obs_arr(3) - obs_arr(1))/inv + 1;
   yrec = abs(obs_arr(4) - obs_arr(2))/inv + 1;
%    obs_arr = [95.1144 -0.9087 105.7139 21.8808]*1e3;
%    nrec = 50;
   
   xl = linspace(obs_arr(1),obs_arr(3),xrec);
   yl = linspace(obs_arr(2),obs_arr(4),yrec);
%    XX = linspace(obs_arr(1),obs_arr(3),nrec);
%    YY = linspace(obs_arr(2),obs_arr(4),nrec);
   
   [XX,YY] = meshgrid(xl,yl);
%    ZZ = 7e3 * ones(size(XX));
   ZZ = depth * 1e3 * ones(size(XX));    % add depth sensitivity test (maybe not used in the output file)
   
   XX = reshape(XX,[],1);
   YY = reshape(YY,[],1);
   ZZ = reshape(ZZ,[],1);
   
   NR = [yrec,xrec];
%    NR = nrec;
   S = layered_wang(xs,ys,zs,strike,dip,rake,U,lp,wp,NR,save_dir,green_dir,count, ...
                   'out_file',out_type,'obs_type',obs_type,'obs_arr',obs_arr);
   
%    fid = fopen([save_dir,'/stress_tensor.out'],'wt');
%    strs_depth = [YY/1e3,XX/1e3,ZZ/1e3,S/1e6];
%    fprintf(fid,'%.4e  %.4e  %.4e  %.4e  %.4e  %.4e  %.4e  %.4e  %.4e\n',strs_depth');
%    fclose(fid);   
end