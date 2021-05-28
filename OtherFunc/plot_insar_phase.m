clc
clear
format long
[track_name,Nsamp]=textread('los_track_L.list','%s %d\n');
Ntrack=length(track_name);

for i=1:Ntrack;
    [x,y,z]=grdread2(['../data/',track_name{i},'.grd']);
    figure;
    h=imagesc(x,y,z);
    set(gca,'YDir','Normal');
    axis equal
    colorbar
    colormap('jet');
    set(h,'alphadata',~isnan(z));
end