function [h0,h1,h2]=plot_insar_sample_new(xinsar,yinsar,zinsar,zout,xx1,xx2,yy1,yy2,varargin)
% plot out InSAR image and its downsampling
% 
% Usage: h=plot_insar_sample(xinsar,yinsar,zinsar,zout,xx1,xx2,yy1,yy2);
% Added fault trace to be plotted
%
% by Kang Wang on 08/27/2015
% by Zeyu Jin on 10/15/2019

format long
set(0,'defaultTextFontName', 'Hevitica');
set(0,'defaultTextFontSize', 12);
Npt=length(zout);
fault_file = [];

if ~isempty(varargin)
   for CC = 1:floor(length(varargin)/2)
      try
         switch lower(varargin{CC*2-1})
             case 'fault'
                fault_file = varargin{CC*2};
         end
      catch
         error('Unrecognized Arguments!');
      end
   end
end

if ~isempty(fault_file)
   fault_trace = load(fault_file);
   lonf = [fault_trace(:,1);fault_trace(:,3)];  
   latf = [fault_trace(:,2);fault_trace(:,4)];
   LS = length(lonf) / 2;
end

XX=[xx1';xx2';xx2';xx1'];
YY=[yy1';yy1';yy2';yy2'];
C=[zout';zout';zout';zout'];

% % too large range
% zmin=min(min(zinsar(~isnan(zinsar)))); 
% zmax=max(max(zinsar(~isnan(zinsar)))); 
% rz=zmax-zmin;
% 
% cmin=zmin+0.5*rz;
% cmax=zmax-0.5*rz;
cmean = nanmean(zinsar(:));
cstd = std(zinsar(:),'omitnan');
cmin = cmean - 8*cstd;
cmax = cmean + 8*cstd;

xmin=min(xinsar);
xmax=max(xinsar);
ymin=min(yinsar);
ymax=max(yinsar);

h0=figure('units','normalized','outerposition',[0 0 1 1]);
set(h0,'renderer','painters');

psv1=[0.05 0.3 0.45 0.55];
h1=subplot('position',psv1);
hold on
h=imagesc(xinsar,yinsar,zinsar);
if ~isempty(fault_file)
   for ii = 1:LS
       line([lonf(ii) lonf(ii+LS)],[latf(ii) latf(ii+LS)],'color','black','linewidth',1.5);
   end
end
set(gca,'YDir','Normal');
set(h,'alphadata',~isnan(zinsar));
axis equal;
hcolor=colorbar('north','Position',[0.1 0.25,0.3,0.01]);
colormap('jet')
caxis([cmin cmax]);
axis([xmin xmax ymin ymax])

psv2=[0.5 0.3 0.45 0.55];
h2=subplot('position',psv2);
hold on
% patch(XX,YY,C);   % plot the scatter points instead
sz = 30;
pxc = (xx1 + xx2) / 2;
pyc = (yy1 + yy2) / 2;
scatter(pxc,pyc,sz,zout,'filled');

% plot the fault segments
if ~isempty(fault_file)
   for ii = 1:LS
       line([lonf(ii) lonf(ii+LS)],[latf(ii) latf(ii+LS)],'color','black','linewidth',1.5);
   end
end

axis equal
hcolor=colorbar('north','Position',[0.55 0.25,0.3,0.01]);
colormap('jet')
caxis([cmin cmax])
axis([xmin xmax ymin ymax])
title(['#pt=',num2str(Npt)]);

end
