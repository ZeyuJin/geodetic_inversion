function create_edcmp_input(slip_model_mat,data_mat,this_dir)

   % d = load('slip_model_for_yuri.mat');
   d = load(slip_model_mat);
   slip_model = d.slip_model;
   xs = slip_model(:,4);
   ys = slip_model(:,5);
   zs = slip_model(:,6);
   lp = slip_model(:,7);
   wp = slip_model(:,8);
   strike = slip_model(:,9);
   dip = slip_model(:,10);
   Np = length(xs);

   % fid = fopen('dir.list');
   % D = textscan(fid,'%s\n');
   % dirlist = D{:,1};  ndata = length(dirlist);
   % fclose(fid);

   %% create input files for unit slip of each patch (strike + dip)
   %% for the edcmp (included in layered_wang.m)

   % for ii = 1:ndata
   % this_dir = dirlist{ii};
   % sp = load([this_dir,'/los_samp3.mat']);
   cmdstr = ['mkdir -p ',this_dir,'/outdata'];  % mkdir for the edcmp output
   system(cmdstr);

   %% each dataset is parallel in CSHELL
   %% function to be compiled
   if contains(this_dir,'GPS')
      sp = load([this_dir,'/gps.mat']);   % GPS data named with "gps.mat"
      xr = sp.data_gps(:,1);
      yr = sp.data_gps(:,2);
   else
      sp = load([this_dir,'/',data_mat]); 
      xr = sp.sampled_insar_data(:,1);
      yr = sp.sampled_insar_data(:,2);
   end

   for jj = 1:Np
       slip = 1;
       rake = 0;
       [uxs,uys,uzs] = layered_wang(xs(jj), ys(jj), zs(jj), strike(jj), dip(jj), rake, slip, lp(jj), wp(jj), xr, yr,this_dir,2*jj-1);
       rake = 90;
       [uxd,uyd,uzd] = layered_wang(xs(jj), ys(jj), zs(jj), strike(jj), dip(jj), rake, slip, lp(jj), wp(jj), xr, yr,this_dir,2*jj);
       U = [uxs,uys,uzs,uxd,uyd,uzd];
       count = 2*jj;
       save([this_dir,'/outdata/3D_disp_',num2str(count),'.mat'],'U');
   end
   % end

end
