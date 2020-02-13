close all
clc
clear

figure; hold on
subplot('Position',[0.07 0.55 0.43 0.4]);
f1 = load('north_dip_misfit.mat');
d1 = f1.misfit_curve;
depth1 = 1:10;
plot(depth1,d1,'r*-','linewidth',2);
xlabel('dipping depth (km)');
xticks([1:10]);
ylabel('RMS reduction (%)');
title('Northern dipping segments');
set(gca,'fontsize',20);
set(gcf,'PaperPositionMode','auto');

subplot('Position',[0.55 0.55 0.43 0.4]);
f2 = load('center_dip_misfit.mat');
d2 = f2.misfit_curve;
depth2 = 1:15;
plot(depth2,d2,'g*-','linewidth',2);
xlabel('dipping depth (km)');
xticks([1:15]);
title('Central dipping segments');
set(gca,'fontsize',20);
set(gcf,'PaperPositionMode','auto');

subplot('Position',[0.25 0.05 0.43 0.4]);
f3 = load('south_dip_misfit.mat');
d3 = f3.misfit_curve;
depth3 = 1:10;
plot(depth3,d3,'b*-','linewidth',2);
xlabel('dipping depth (km)');
xticks([1:10]);
ylabel('RMS reduction (%)');
title('Southern dipping segments');
set(gca,'fontsize',20);
set(gcf,'PaperPositionMode','auto');
