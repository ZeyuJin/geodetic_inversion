function [G_raw,G,bdata_raw,bdata] = build_green_function(slip_model,sampled_data_file,option,ramp_choice,model_type,varargin)
% build the Green's function and data column based on slip model and
% resampled InSAR data, and assume uniform weight factors
% according to Simons et al., BSSA 2002
% model_type could be 'okada' or 'layered' (EDGRN/EDCMP)
   
   data = load(sampled_data_file);
   if ~isempty(varargin)
      beta = varargin{1};
   else
      beta = 0.3;   % relative weighting between GPS and InSAR as a default
   end
   ramp_choice = lower(ramp_choice);
   
   % for a layered model
   [filepath,~,~] = fileparts(sampled_data_file);
   
   switch option
       case {'insar','AZO'}
           sampled_data = data.sampled_insar_data;
           h1 = size(sampled_data,1);
           
%            rms_insar = double(data.rms_out); 
           rms_insar = ones(h1,1);      % uniform weighting for InSAR/AZO data
           w_insar = calc_weight_insar_error(rms_insar);
           
           if strcmp(option,'insar')
               if strcmp(model_type,'okada')
                   G_raw = calc_green_insar_okada(slip_model,sampled_data);
               else
                   GG = load([filepath,'/layered_green.mat']);        % from EDCMP output
%                    GG = load([filepath,'/modelB_layer.mat']);
                   G_raw = GG.G_raw;
               end
           else
               if strcmp(model_type,'okada')
                   G_raw = calc_green_AZO_okada(slip_model,sampled_data);
               else
                   GG = load([filepath,'/layered_green.mat']);
%                    GG = load([filepath,'/modelB_layer.mat']);
                   G_raw = GG.G_raw;
               end
           end
           bdata_raw = sampled_data(:,3); 
           
           if strcmp(ramp_choice,'bi_ramp')        % assume bilinear ramp 
               xsar = sampled_data(:,1) / 1000;        % in km
               ysar = sampled_data(:,2) / 1000;        % in km
               dem_out = data.dem_out / 1000;          % in km               
               rmp = [xsar,ysar,dem_out,ones(size(xsar))];
           elseif strcmp(ramp_choice,'qu_ramp_7')
               xsar = sampled_data(:,1) / 1000;        % in km
               ysar = sampled_data(:,2) / 1000;        % in km
               dem_out = data.dem_out / 1000;          % in km  
               rmp = [xsar.^2,ysar.^2,xsar.*ysar,xsar,ysar,dem_out,ones(size(xsar))];
           elseif strcmp(ramp_choice,'qu_ramp_5')
               xsar = sampled_data(:,1) / 1000;        % in km
               ysar = sampled_data(:,2) / 1000;        % in km
               dem_out = data.dem_out / 1000;          % in km  
               rmp = [xsar.*ysar,xsar,ysar,dem_out,ones(size(xsar))];
           else
               rmp = [];
           end
           
           G_raw = double([G_raw,rmp]);
           G = w_insar * G_raw;
           bdata = w_insar * bdata_raw;            % in cm
           
       case 'cgps'
           data_gps = data.data_gps;
           Ngps = 3 * size(data_gps,1);            % three components
           
           if strcmp(model_type,'okada')
               G_raw = calc_green_gps_3d_okada(slip_model,data_gps);
           else
               GG = load([filepath,'/layered_green.mat']);
%                GG = load([filepath,'/modelB_layer.mat']);
               G_raw = GG.G_raw;
           end
               
           bdata_raw = double([data_gps(:,3);data_gps(:,4);data_gps(:,5)]);          
           sig_gps = double([data_gps(:,6);data_gps(:,7);data_gps(:,8)]);
           w_gps = calc_weight_gps_error(sig_gps);
           w_gps = w_gps .* beta;   
           
           if strcmp(ramp_choice,'bi_ramp')
               rmp = zeros(Ngps,4);                % assume bilinear ramp
           elseif strcmp(ramp_choice,'qu_ramp_7')
               rmp = zeros(Ngps,7);
           elseif strcmp(ramp_choice,'qu_ramp_5')
               rmp = zeros(Ngps,5);
           else
               rmp = [];
           end
           
           G_raw = double([G_raw,rmp]);
           G = w_gps * G_raw;
           bdata = w_gps * bdata_raw;
           
       case 'camp_gps'   % only use horizontal displacements
           data_gps = data.data_gps;
           Ngps = 2 * size(data_gps,1);
           
           if strcmp(model_type,'okada')
               G_raw = calc_green_gps_2d_okada(slip_model,data_gps);
           else
               GG = load([filepath,'/layered_green.mat']);
%                GG = load([filepath,'/modelB_layer.mat']);
               G_raw = GG.G_raw;
           end
           
           bdata_raw = double([data_gps(:,3);data_gps(:,4)]);
           sig_gps = double([data_gps(:,5);data_gps(:,6)]);
           w_gps = calc_weight_gps_error(sig_gps);
           w_gps = w_gps .* beta;
           
           if strcmp(ramp_choice,'bi_ramp')
               rmp = zeros(Ngps,4);                % assume bilinear ramp
           elseif strcmp(ramp_choice,'qu_ramp_7')
               rmp = zeros(Ngps,7);
           elseif strcmp(ramp_choice,'qu_ramp_5')
               rmp = zeros(Ngps,5);
           else
               rmp = [];
           end
           
           G_raw = double([G_raw,rmp]);
           G = w_gps * G_raw;
           bdata = w_gps * bdata_raw;
   end 
end