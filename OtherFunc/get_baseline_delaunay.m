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

figure;
plot(t_sar,Bp,'ro','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','r')
text(t_sar,Bp+3,date_sar,'HorizontalAlignment','center')
axis([xmin xmax ymin ymax])
hold on
%plot([0 0],[ymin+5 ymax-5],'r--','LineWidth',2);
xlabel('Time (days since the earthquake)');
ylabel('Perpendicular Baseline (m)');
dt=delaunayTriangulation(t_sar,Bp);
TRI=delaunay(t_sar,Bp);
ntri=length(TRI(:,1));
nside=ntri*3;  %total number of sides, including repeating ones


indx_master=[];
indx_slave=[];
for i=1:ntri;
    this_element=TRI(i,:);
    indx_a=this_element(1);
    indx_b=this_element(2);
    indx_c=this_element(3);
    indx_element=[indx_a,indx_b,indx_c];
    
        t1=t_sar(indx_a);
        t2=t_sar(indx_b);
        if (t1<t2);
           indx_master_tmp=indx_a;
           indx_slave_tmp=indx_b;
        else
            indx_master_tmp=indx_b;
            indx_slave_tmp=indx_a;
        end
        
        
        indx_master=[indx_master;indx_master_tmp];
        indx_slave=[indx_slave;indx_slave_tmp];
        
        t1=t_sar(indx_b);
        t2=t_sar(indx_c);
        if (t1<t2);
           indx_master_tmp=indx_b;
           indx_slave_tmp=indx_c;
        else
            indx_master_tmp=indx_c;
            indx_slave_tmp=indx_b;
        end
        
        indx_master=[indx_master;indx_master_tmp];
        indx_slave=[indx_slave;indx_slave_tmp];
        
        t1=t_sar(indx_c);
        t2=t_sar(indx_a);
        if (t1<t2);
           indx_master_tmp=indx_c;
           indx_slave_tmp=indx_a;
        else
            indx_master_tmp=indx_a;
            indx_slave_tmp=indx_c;
        end
        indx_master=[indx_master;indx_master_tmp];
        indx_slave=[indx_slave;indx_slave_tmp];
        
    

end

indx_intf=[indx_master,indx_slave];
indx_intf_unique=unique(indx_intf,'rows');
Ninsar=length(indx_intf_unique(:,1));

f_insar=fopen('intf_delaunay.in','w');

for i=1:Ninsar;
    indx_this_intf=indx_intf_unique(i,:);
    indx_first=indx_this_intf(1);
    indx_second=indx_this_intf(2);
    t_this_intf=t_sar([indx_first,indx_second]);
    Bp_this_intf=Bp([indx_first,indx_second]);
    scene1=scene_id{indx_first};
    scene2=scene_id{indx_second};
    
    scene_intf=[scene1,':',scene2];
    
    plot(t_this_intf,Bp_this_intf,'b-');
    fprintf(f_insar,'%s\n',scene_intf);
end

fclose(f_insar);

