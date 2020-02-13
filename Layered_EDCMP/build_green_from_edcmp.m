function build_green_from_edcmp(slip_model_mat,data_list)

    fid = fopen(data_list);
    D = textscan(fid,'%s %s\n');
    dirlist = D{:,1};  TYPE = D{:,2};
    ndata = length(dirlist);
    fclose(fid);
   
    d = load(slip_model_mat);
    slip_model = d.slip_model;
    NP = size(slip_model,1);
    
%    fid = fopen('all_data');
%    C = textscan(fid,'%s  %s\n');
%    data_path = C{:,1};    TYPE = C{:,2};
%    fclose(fid);
    
    for ii = 1:ndata
        this_track = dirlist{ii};
        data_type = TYPE{ii};
        if strcmp(data_type,'insar') || strcmp(data_type,'AZO')
            data = load([this_track,'/los_samp_detrend_mask.mat']);   % mask out the near field data
            ve = data.sampled_insar_data(:,4);
            vn = data.sampled_insar_data(:,5);
            vz = data.sampled_insar_data(:,6);
            G_raw = zeros(length(ve),2*NP);
            for jj = 1:NP
                count = jj*2;
                load([this_track,'/outdata/3D_disp_',num2str(count),'.mat']);
%                 U = U .* 100;  % convert to cm
                if strcmp(data_type,'insar')                    
                    U_los1 = U(:,1).*ve + U(:,2).*vn + U(:,3).*vz;
                    U_los2 = U(:,4).*ve + U(:,5).*vn + U(:,6).*vz;
                    G_raw(:,jj) = U_los1;
                    G_raw(:,jj+NP) = U_los2;
                else
                    theta_az = -atan2d(vn,ve) - 180;
                    U_az1=U(:,1).*sind(theta_az) + U(:,2).*cosd(theta_az);
                    U_az2=U(:,4).*sind(theta_az) + U(:,5).*cosd(theta_az);
                    G_raw(:,jj) = U_az1;
                    G_raw(:,jj+NP) = U_az2;                    
                end
            end
        elseif strcmp(data_type,'cgps')    % continuous GPS data
            data = load([this_track,'/gps.mat']);
            Ngps = length(data.data_gps(:,1));
            G_raw = zeros(3*Ngps,2*NP);
            for jj = 1:NP
                count = jj*2;   % fix the bug
                load([this_track,'/outdata/3D_disp_',num2str(count),'.mat']);
%                 U = U .* 100;  % convert to cm
                U_g1 = [U(:,1);U(:,2);U(:,3)];
                U_g2 = [U(:,4);U(:,5);U(:,6)];
                G_raw(:,jj) = U_g1;
                G_raw(:,jj+NP) = U_g2;
            end
        elseif strcmp(data_type,'sgps')   % survey GPS data
            data = load([this_track,'/gps.mat']);
            Ngps = length(data.data_gps(:,1));
            G_raw = zeros(2*Ngps,2*NP);
            for jj = 1:NP
                count = jj*2;   % fix the bug
                load([this_track,'/outdata/3D_disp_',num2str(count),'.mat']);
                U_g1 = [U(:,1);U(:,2)];    % only consider the horizontal components
                U_g2 = [U(:,4);U(:,5)];
                G_raw(:,jj) = U_g1;
                G_raw(:,jj+NP) = U_g2;
            end           
        else
            error('Data type is wrong and exit!');
        end
        save([this_track,'/layered_green.mat'],'G_raw');
     end
end
