function show_slip_fem(slip_model);
%function: display the the slip vectors on the node
% by Kang Wang on 08/07/2014
% Last Updated by Kang Wang on 02/11/2016
%
%  Usage: show_slip_fem(data);
%  where:
%     id_node     = data(:,1);
%     xnode       = data(:,2);
%     ynode       = data(:,3);
%     znode       = data(:,4);
%     strike_node = data(:,5);
%     strike_x    = data(:,6);
%     dip_node    = data(:,7);
%     slip1       = data(:,8);
%     slip2=      = data(:,9);


%set(0,'defaultAxesFontName', 'AvantGarde')
%set(0,'defaultAxesFontSize', 15)
d2r=pi/180.0;
%data=load(slip_model);
data=slip_model;
id_node=data(:,1);  %id of the node
xnode=data(:,2);
ynode=data(:,3);
znode=data(:,4);
strike_node=data(:,5);
strike_x=data(:,6);
dip_node=data(:,7);
slip1=data(:,8);
slip2=data(:,9);

slip=sqrt(slip1.^2+slip2.^2);

% strike_x: strike of the x-axis of the FEM MODEL
strike_x=mean(strike_x);
dtheta_flt=(strike_node-strike_x)*d2r;
st_H_flt=slip2.*cos(dip_node*d2r);
sx_flt=slip1.*cos(dtheta_flt)+st_H_flt.*sin(dtheta_flt);
sz_flt=slip2.*sin(dip_node*d2r);
sy_flt=-slip2.*sin(dtheta_flt)+st_H_flt.*cos(dtheta_flt);

if (sqrt(sum(xnode.^2))<10);
hs=scatter(ynode/1000,znode/1000,25,slip,'filled');
alpha(hs,0.5)
%colorbar
axis equal
hold on
grid on
quiver(ynode/1000,znode/1000,sy_flt,sz_flt,...
 'color',[0 0 0],'LineWidth',0.5);
%xlabel(['Y-FEM (km) strike=',num2str(strike_x+90),' deg']);
%ylabel(['Depth (km)']);
elseif (sqrt(sum(ynode.^2))<10);
scatter(xnode/1000,znode/1000,25,slip,'filled');
%colorbar
axis equal
hold on
grid on
quiver(xnode/1000,znode/1000,sx_flt,sz_flt,...
 'color',[0 0 0],'LineWidth',0.);
%xlabel(['X-FEM (km) strike=',num2str(strike_x),' deg']);
%ylabel(['Depth (km)']);
elseif (sqrt(sum(znode.^2))<10);
scatter(xnode/1000,ynode/1000,25,slip,'filled');
%colorbar
axis equal
hold on
grid on
quiver(xnode/1000,ynode/1000,sx_flt,sy_flt,...
 'color',[0 0 0],'LineWidth',0.);
%xlabel(['X-FEM (km) strike=',num2str(strike_x),' deg']);
%ylabel(['Y-FEM (km)']);
else
scatter3(xnode/1000,ynode/1000,znode/1000,25,slip,'filled');
%colorbar
axis equal
hold on
grid on
quiver3(xnode/1000,ynode/1000,znode/1000,sx_flt,sy_flt,1.1*sz_flt,...
       'color',[0 0 0],'LineWidth',0.5);
%xlabel(['X-FEM(km) strike=',num2str(strike_x),' deg']);
%ylabel(['Y-FEM (km)']);
end
axis tight
