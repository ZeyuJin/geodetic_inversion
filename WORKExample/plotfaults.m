data=importdata('all_faults');
[x0,y0]=ll2xy(37,37,37);
% newdata=zeros(1,4);
% for i=1:12
%     [x1,y1]=ll2xy(data(i,1),data(i,2),37);
%     [x2,y2]=ll2xy(data(i,3),data(i,4),37);
%     x1=(x1-x0)/1000;
%     y1=(y1-y0) / 1000;
%     x2= (x2-x0)/1000;
%     y2 = (y2-y0) / 1000;
%     newdata=[newdata;[x1,y1,x2,y2]];
% end
% newdata(1,:)=[];
for i=1:12
    plot([data(i,1) data(i,3)],[data(i,2) data(i,4)])
    text((data(i,1)+data(i,3))/2,(data(i,2)+data(i,4))/2,num2str(i))
    hold on
end
    
%plot([data(:,1) data(:,3)],[data(:,2) data(:,4)]);

%possible dipping faults: #8 #6 