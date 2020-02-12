function detrend_range_offsets(filepath,losfile,rngfile,DEMfile,ramp_type,varargin)
% detrend the range offsets using LOS data
% Subtract the range offsets from LOS displacements, and then 
% detrend the residual left
% save with the grid name: los_clean_detrend.grd
% the los,rng,DEM grids should be in the same size 

   this_track = filepath;  
   lon_eq = -117.5;
   lat_eq = 35.5;
   res_max = 50;
   
   %% read varargin values and assembly
   if ~isempty(varargin)
       for CC = 1:floor(length(varargin)/2)
           try
               switch lower(varargin{CC*2-1})
                   case 'lonc'
                       lon_eq = varargin{CC*2};
                   case 'latc'
                       lat_eq = varargin{CC*2};
                   case 'misfit_range'
                       res_max = varargin{CC*2};
               end
           catch
               error('Unrecognized Keyword\n');
           end
       end
   end
   
   [xo,yo] = utm2ll(lon_eq,lat_eq,0,1);  
   [~,~,demin] = grdread2([this_track,'/',DEMfile]);       % usually dem_samp.grd
   [~,~,losin] = grdread2([this_track,'/',losfile]);       % usually los_cut.grd
   [lon,lat,rngin] = grdread2([this_track,'/',rngfile]);   % usually unwrap_sample_clean.grd
   res = rngin - losin;
   
   [mlon,mlat] = meshgrid(lon,lat);
   [xsar,ysar] = utm2ll(mlon(:),mlat(:),0,1);
   xsar = (xsar - xo) ./ 1000;   
   ysar = (ysar - yo) ./ 1000;
   demin = demin ./ 1000;
   
   xsar = reshape(xsar,size(mlon));
   ysar = reshape(ysar,size(mlat));

   if strcmp(ramp_type,'bi_ramp')
      pfit = fit_ramp_topo(double(res),xsar,ysar,demin);
      ramp = pfit(1).*xsar + pfit(2).*ysar + pfit(3).*demin + pfit(4);
      
   elseif strcmp(ramp_type,'qu_ramp_7')
      pfit = fit_ramp_topo_quadratic(double(res),xsar,ysar,demin);
      ramp = pfit(1).*xsar.^2 + pfit(2).*ysar.^2 + pfit(3).*xsar.*ysar + pfit(4).*xsar + ...
             pfit(5).*ysar + pfit(6).*demin + pfit(7);
         
   elseif strcmp(ramp_type,'qu_ramp_5')
      pfit = fit_ramp_topo_hyper(double(res),xsar,ysar,demin);
      ramp = pfit(1).*xsar.*ysar + pfit(2).*xsar + pfit(3).*ysar + pfit(4).*demin + pfit(5);
      
   elseif strcmp(ramp_type,'sine')
      ramp = fit_sinusoid(double(res),xsar,ysar);
      
   else
      ramp = zeros(size(rngin));     % to be improved
   end

   
   rng_detrend = rngin - ramp;
   res_detrend = res - ramp;
   figure;
   subplot('Position',[0.03 0.55 0.45 0.4]);
   pcolor(xsar,ysar,rngin);
   shading flat
   colormap jet
   colorbar
   title('Original offsets (cm)');
   set(gca,'Fontsize',20);
       
   subplot('Position',[0.53 0.55 0.45 0.4]);
   pcolor(xsar,ysar,rng_detrend);
   shading flat
   colormap jet
   colorbar
   title('De-trended offsets (cm)');
   set(gca,'Fontsize',20);
   
   subplot('Position',[0.03 0.05 0.45 0.4]);
   pcolor(xsar,ysar,res);
   shading flat
   colormap jet
   colorbar
   title('Residual of offset and LOS (cm)');
   set(gca,'Fontsize',20);
   caxis([0 res_max]);
       
   subplot('Position',[0.53 0.05 0.45 0.4]);
   pcolor(xsar,ysar,res_detrend);
   shading flat
   colormap jet
   colorbar
   title('Fitting quadratic ramp (cm)');
   set(gca,'Fontsize',20);
   caxis([-20 20]);
   
   set(gcf,'PaperPositionMode','auto');
   grdwrite2(lon,lat,rng_detrend,[this_track,'/los_clean_detrend.grd']);

end
