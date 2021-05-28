function [xout,yout]=get_polygon(fig_id,N);

xout=zeros(N+1,1);
yout=zeros(N+1,1);

figure(fig_id);
hold on
for i=1:N;
  [xpt,ypt]=ginput(1);
   xout(i)=xpt;
   yout(i)=ypt;
   plot(xpt,ypt,'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','r');
   title(['Nclick: ',num2str(i),'/',num2str(N)]);
end
xout(N+1)=xout(1);
yout(N+1)=yout(1);

plot(xout,yout,'b-','LineWidth',1.5)
