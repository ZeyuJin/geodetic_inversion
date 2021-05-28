function fix_snaphu_ambiguity(dir_data,grd_list);
%
% Estimate the phase ambiguity between adjcent frames or subswaths 
% due to snaphu
%
% Usage: fix_snaphu_ambiguity(dir_data,grd_list);
%
% 
% Examples of contents in grd_list;
%  unwrap_IW1_mask.grd
%  unwrap_IW2_mask.grd
%  unwrap_IW3_mask.grd  
%
% This function should pop up a window displaying the consecutive
% images and and wait for clicking points to define region of overlapping.
% Currently, 6 points are required to define the polygon of overlapping
% If more points needed, change Nclick in the code accordingly. 
% 
%
% plase make sure data of consecutive lines have overlapping
% The output of the function of .grd files by adding integer of 
% 2*pi to the orginal .grd files and be renamed to *_new.grd
%
% After running this, please go to dir_data and run 
%   merge_grd.sh merge.list
% to generate the merged .grd file.
%
% by Kang Wang in Nov. 2015
 
format long
snaphu_list=textread(grd_list,'%s');
Nsnaphu=length(snaphu_list);

for i=1:Nsnaphu-1;
  name1=snaphu_list{i};
  name2=snaphu_list{i+1};
  root1=name1(1:length(name1)-4);
  root2=name2(1:length(name2)-4);
   
  if (i==1);
    grd1_in=[dir_data,root1,'.grd'];
    grd2_in=[dir_data,root2,'.grd'];
  else
    grd1_in=[dir_data,root1,'_new.grd'];
    grd2_in=[dir_data,root2,'.grd'];
  end

  [x1,y1,z1]=grdread2(grd1_in);
  [x2,y2,z2]=grdread2(grd2_in);
  [X1,Y1]=meshgrid(x1,y1);
  [X2,Y2]=meshgrid(x2,y2);

  indx_good1=~isnan(z1);
  indx_good2=~isnan(z2);

  xgood1=X1(indx_good1);
  ygood1=Y1(indx_good1);
  zgood1=z1(indx_good1);

  xgood2=X2(indx_good2);
  ygood2=Y2(indx_good2);
  zgood2=z2(indx_good2);

  pos1=[xgood1,ygood1];
  pos2=[xgood2,ygood2];

  x1min=min(xgood1);
  x1max=max(xgood1);
  y1min=min(ygood1);
  y1max=max(ygood1);

  x2min=min(xgood2);
  x2max=max(xgood2);
  y2min=min(ygood2);
  y2max=max(ygood2);
  xmin=min([x1min x2min]);
  xmax=max([x1max x2max]);
  ymin=min([y1min y2min]);
  ymax=max([y1max y2max]);

  figure('units','normalized','outerposition',[0 0 1 1])
  h1=imagesc(x1,y1,z1);
  set(h1,'alphadata',~isnan(z1));
  hold on
  h2=imagesc(x2,y2,z2);
  set(h2,'alphadata',~isnan(z2));
  set(gca,'YDir','Normal')
  axis equal
  axis([xmin xmax ymin ymax])
  colorbar
  colormap('jet')
  title([root1,'---',root2],'interpreter','none');
  Nclick=6;
  lon_pt=zeros(Nclick,1);
  lat_pt=zeros(Nclick,1);
  for k=1:Nclick;
     [lon_pt(k),lat_pt(k)]=ginput(1);
     plot(lon_pt(k),lat_pt(k),'ko','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','g');
  end

  lon_plot=[lon_pt;lon_pt(1)];
  lat_plot=[lat_pt;lat_pt(1)];
  plot(lon_plot,lat_plot,'k-','LineWidth',2);
  
  in1=inpolygon(xgood1,ygood1,lon_plot,lat_plot);
  xfind1=xgood1(in1);
  yfind1=ygood1(in1);
  zfind1=zgood1(in1);
  
  in2=inpolygon(xgood2,ygood2,lon_plot,lat_plot);
  xfind2=xgood2(in2);
  yfind2=ygood2(in2);
  zfind2=zgood2(in2);
  
  
  idx=knnsearch([xfind2,yfind2],[xfind1,yfind1]);
  xxfind1=xfind1;
  yyfind1=yfind1;
  zzfind1=zfind1;
  
  xxfind2=xfind2(idx);
  yyfind2=yfind2(idx);
  zzfind2=zfind2(idx);
  
  dz=zzfind2-zzfind1;
%  N2pi=round(mean(dz)/2/pi);
  N2pi=mean(round(dz/2/pi));
  
  z1_new=z1;
  z2_new=z2-N2pi*2*pi;




  
  grd1_out=[dir_data,root1,'_new.grd'];
  grd2_out=[dir_data,root2,'_new.grd'];
  grdwrite2(x1,y1,z1_new,grd1_out);
  grdwrite2(x2,y2,z2_new,grd2_out);
  
end

f_blend=fopen([dir_data,'merge.list'],'w');
figure;
hold on
for i=1:Nsnaphu;
  name_out=snaphu_list{i};
  root_out=name_out(1:length(name_out)-4);

  grd_out=[dir_data,root_out,'_new.grd']; 
  fprintf(f_blend,'%s\n',[root_out,'_new.grd']);
  [xout,yout,zout]=grdread2(grd_out);

  h=imagesc(xout,yout,zout);
  set(h,'alphadata',~isnan(zout));
end
set(gca,'YDir','Normal')
axis equal
colorbar 
colormap('jet');

fclose(f_blend);
