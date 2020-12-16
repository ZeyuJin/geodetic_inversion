function los_model = slip2disp_edcmp(this_track,xin,yin,losl,zel,znl,zul,slip_model,data_type)
% this script will compute the layered prediction using EDCMP
% based on the full resolution data
% not the sub-sampled scatter points

   % save the binary irregular observation points
   temp = [yin(:),xin(:)]';   % x is north in EDCMP
   fid = fopen([this_track,'/full_obs.rec'],'w');
   fwrite(fid,temp,'real*4');
   fclose(fid);
   NR = length(yin(:));   % number of observation points
   
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
   cmdstr = ['mkdir -p ',this_track,'/outdata'];  % mkdir for the edcmp output
   system(cmdstr);      
   count = '';
   U = layered_wang(xs,ys,zs,strike,dip,rake,U,lp,wp,NR,this_track,count);
%    out = load([this_track,'/outdata/data_only']); 
%    ux=out(:,4);
%    uy=out(:,3);
%    uz=-out(:,5);
   
   % convert to matrix form and LOS/AZO 
   ux = reshape(U(:,1),size(losl));
   uy = reshape(U(:,2),size(losl));
   uz = reshape(U(:,3),size(losl));
   
   if strcmp(data_type,'insar')
       los_model = ux.*zel + uy.*znl + uz.*zul;
   elseif strcmp(data_type,'AZO')
       theta_az = -atan2d(znl,zel) - 180;
       los_model = ux.*sind(theta_az) + uy.*cosd(theta_az);
   elseif strcmp(data_type,'3d')
       los_model = [ux,uy,uz];
   else
       error('There is something wrong with data type!');
   end
   los_model = los_model .* 100;   % convert to cm finally
   
end