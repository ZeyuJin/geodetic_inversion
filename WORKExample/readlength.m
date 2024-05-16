clear
clc
faults=importdata('all_faults');
lon_c=37;
newcoor=zeros(1,5);
for i=1:5
    [xo1,yo1]=ll2xy(faults(i,1),faults(i,2),lon_c);
    [xo2,yo2]=ll2xy(faults(i,3),faults(i,4),lon_c);
    len=sqrt((xo1-xo2)^2+(yo1-yo2)^2);
    newcoor=[newcoor;[xo1,yo1,xo2,yo2,len]];
end
newcoor(1,:)=[];