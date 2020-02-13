function gps_model2ascii(gps_model_mat,Ngps,model_type,data_type)
    load(gps_model_mat);
    modelx = gps_model(1:Ngps);
    modely = gps_model(Ngps+1:2*Ngps);
    model = [modelx,modely];
    
    dd = load(['GPS/',data_type,'_gps_data_error.txt']);
    lon = dd(:,1);
    lat = dd(:,2);
    error = zeros(length(lon),3);
    mat = [lon,lat,model,error]';
    
    if strcmp(model_type,'okada')
        save_name = [data_type,'_gps_model.txt'];
    else
        save_name = [data_type,'_gps_model_layer.txt'];
    end
    fid = fopen(['GPS/',save_name],'w');
    fprintf(fid,'%f  %f  %f  %f  %f  %f  %f\n',mat);
    fclose(fid); 
end