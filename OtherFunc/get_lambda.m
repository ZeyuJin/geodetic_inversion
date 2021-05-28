function lambda_find=get_lambda(data_lambda);
%
% plot out the L_curve and let the user select the point on the plot 
% 
% Usage: lambda_find=get_lambda(data_lambda);
%
% content of data_lambda:
% 
%   data_lambda=[lambda,chi_data,chi_model];
%   
% by Kang Wang in Feb. 2016

lambda_in=data_lambda(:,1);
chi_data=data_lambda(:,2);
chi_model=data_lambda(:,3);

[lambda_sort,indx_sort]=sort(lambda_in,'descend');
chi_data_sort=chi_data(indx_sort);
chi_model_sort=chi_model(indx_sort);

xx=chi_model_sort;
yy=chi_data_sort;

xx_norm=xx/max(xx);
yy_norm=yy/max(yy);

h=figure;
set(h,'render','painters');
plot(xx,yy,'k-','LineWidth',1);
hold on
plot(xx,yy,'ro','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k');
xlabel('Roughness (cm^2/km^2)');
ylabel('\chi^2');
hold on
[xpt,ypt]=ginput(1);

xpt_norm=xpt/max(xx);
ypt_norm=ypt/max(yy);

idx=knnsearch([xx_norm,yy_norm],[xpt_norm,ypt_norm]);
xfind=xx(idx);
yfind=yy(idx);
lambda_find=lambda_sort(idx);

f_lambda=fopen('lambda_find.dat','w');
 fprintf(f_lambda,'%16.6e\n',lambda_find);
fclose(f_lambda);

rx=max(xx)-min(xx);
ry=max(yy)-min(yy);
plot(xfind,yfind,'ro','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r');

str=['$\lambda=$',num2str(lambda_find)];
text(xfind+0.05*rx,yfind+0.05*ry,str,'FontSize',15,'Interpreter','latex');
