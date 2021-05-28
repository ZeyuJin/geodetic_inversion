function gene_slip2cubit_PRM(PRMFILE);
% generate a PRM template for slip2cubit
%
% usage:  gene_slip2cubit_PRM('slip2cubit.PRM');
%
% by Kang Wang in Feb. 2016

f_prm=fopen(PRMFILE,'w');
fprintf(f_prm,'%s\n',['###### PRMFILE FOR SLIP2CUBIT #########']);
fprintf(f_prm,'%s %12.4f\n',['LON_ORIGIN   '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['LAT_ORIGIN   '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['LON_CENTER   '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['STRIKE_X   '],9999.9999);

fprintf(f_prm,'%s\n',['####################################### ']);
fprintf(f_prm,'%s %s\n',['MODEL_ID    '],['AAAAAA']);
fprintf(f_prm,'%s %s\n',['CUBIT_DIR    '],['AAAAAA']);
fprintf(f_prm,'%s %s\n',['FAULT_TRACE_LL    '],['AAAAAA']);
fprintf(f_prm,'%s %s\n',['FAULT_XYZ    '],['AAAAAA']);
fprintf(f_prm,'%s %s\n',['DEM_FILE    '],['AAAAAA']);
fprintf(f_prm,'%s %s\n',['INCLUDE_TOPO    '],['AAAAAA']);
fprintf(f_prm,'%s %s\n',['SLIP_MODEL_OKADA    '],['AAAAAA']);
fprintf(f_prm,'%s\n',['#######################################']);

fprintf(f_prm,'%s %12.4f\n',['XMIN    '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['XMAX    '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['YMIN    '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['YMAX    '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['ZMIN    '],9999.9999);

fprintf(f_prm,'%s\n',['#######################################']);


fprintf(f_prm,'%s\n',['########## TOPOGRAPHY RESOLUTION #######']);
fprintf(f_prm,'%s %12.4f\n',['XMIN_FINE    '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['XMAX_FINE    '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['YMIN_FINE    '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['YMAX_FINE    '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['DX_FINE    '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['DY_FINE    '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['DX_COARSE    '],9999.9999);
fprintf(f_prm,'%s %12.4f\n',['DY_COARSE    '],9999.9999);
disp(['**************** IMPORTANT ***********'])
disp([' This is just a template, please change corresponding values!  '])

fclose(f_prm);

