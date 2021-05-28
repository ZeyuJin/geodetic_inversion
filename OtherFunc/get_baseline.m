clc
clear
set(0,'defaultAxesFontName', 'Helvetica')
set(0,'defaultAxesFontSize', 12)


baseline_table ='baseline_F3_clean.dat';
date_eq='20150425';
day_eq=datenum(date_eq,'yyyymmdd');

[scene_id,sar_time,sar_day,B_para,Bp]=textread(baseline_table,'%s %f %d %f %f\n');
Nsar=length(scene_id);

t_sar=[];
for i=1:Nsar;
   id_this_scene=scene_id{i};
   date_this_scene=char(id_this_scene(4:11));
   date_sar{i}=date_this_scene;
   day_this_scene=datenum(date_this_scene,'yyyymmdd');
   iday_this_scene=day_this_scene-day_eq;
   t_sar=[t_sar;iday_this_scene];
   
end

Rx=max(t_sar)-min(t_sar);
xmin=min(t_sar)-0.05*Rx;
xmax=max(t_sar)+0.05*Rx;

Ry=max(Bp)-min(Bp);
ymin=min(Bp)-0.05*Ry;
ymax=max(Bp)+0.05*Ry;


dt=delaunayTriangulation(t_sar,Bp);
TRI=delaunay(t_sar,Bp);
ntri=length(TRI(:,1));
nside=ntri*3;  %total number of sides, including repeating ones

figure;
plot(t_sar,Bp,'ro','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','r')
text(t_sar,Bp+3,date_sar,'HorizontalAlignment','center')
axis([xmin xmax ymin ymax])
hold on
%plot([0 0],[ymin+5 ymax-5],'r--','LineWidth',2);
xlabel('Time (days since the earthquake)');
ylabel('Perpendicular Baseline (m)');


DBp_max=100;  %maximum baseline allowed
Dt_min=0;     %minmum time interval allowed
Dt_max=60;    %maximum time interval allowded (days) 

Ninsar=0;

f_insar=fopen('intf.in','w');
for i=1:Nsar;
    t_this_sar=t_sar(i);
    Bp_this_sar=Bp(i);
    Dt=t_sar-t_this_sar;
    DBp=abs(Bp-Bp_this_sar);
    
    indx_good=find(Dt>Dt_min & Dt<Dt_max &DBp<DBp_max);
    
    Nfind=length(indx_good);
    if (Nfind>0);
        for n=1:Nfind;
            Ninsar=Ninsar+1;
            t_find=t_sar(indx_good(n));
            Bp_find=Bp(indx_good(n));
            fprintf(f_insar,'%s\n',[scene_id{i},':',scene_id{indx_good(n)}]);
            plot([t_this_sar,t_find],[Bp_this_sar,Bp_find],'b');
        end
    end
end

Ninsar
fclose(f_insar);