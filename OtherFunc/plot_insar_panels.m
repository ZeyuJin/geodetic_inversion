function plot_insar_panels(x,y,z,title_str,Nx,Ny,cmin,cmax,dir_figure_out,figure_id);
% plot out panels for multiple InSAR observations
%
% Usage: plot_insar_panels(x,y,z,title_str,Nx,Ny,cmin,cmax,dir_figure_out,figure_id);
% x,y,z,title_str ----- cells 
% cmin/cmax defines the color range
% 
% figure_id is the output name 
% dir_figure_out is the output directory (e.g. '../figure_test/')/
%
% by Kang Wang on 01/12/2017


Nintf=length(z);
disp(['Number of total images to display: ',num2str(Nintf)]);
if (exist(dir_figure_out)~=7);
  mkdir(dir_figure_out);
end 

disp(['Saving Figures to ',dir_figure_out]);
Npage=Nx*Ny;
dx=0.8/Nx;
dy=0.8/Ny;
nsubplot=0;
for i=1:Ny;
 for j=1:Nx;
   nsubplot=nsubplot+1;
   psv{nsubplot}=[0.1+(j-1)*dx 0.1+(Ny-i)*dy dx-0.01 dy-0.01];
 end 
end

k=0;
nplot=0;

nfigure=0;
for i=1:Nintf;
  if (mod(i,Npage)==0);
    k=k+1;
    hf= figure('units','normalized','outerposition',[0 0 1 1]);
    set(hf,'Visible','off')
    nfigure=nfigure+1;
    for nn=1:Npage;
       mm=(k-1)*Npage+nn;
    if (mm <=Nintf)
    this_x=x{mm};
    this_y=y{mm};
    this_z=z{mm};
    xmax=max(this_x);
    xmin=min(this_y);

    ymax=max(this_y);
    ymin=min(this_y);

    Rx=xmax-xmin;
    Ry=ymax-ymin;
     subplot('position',psv{nn});
     hi=imagesc(this_x,this_y,this_z);
     colormap(jet);
     caxis([cmin cmax])
     set(gca,'XAxisLocation','top')
    set(hi,'alphadata',~isnan(this_z));
%    axis equal
%    axis tight

    if (nn==1);
     hcolor=colorbar('south','Position',[0.1+(Nx-1)*dx 0.915,dx-0.01,0.01]);
    end
    if (nn~=1);
      set(gca,'YTick',[])
      set(gca,'XTick',[])
    end
     set(gca,'YDir','normal')
      text(xmin+0.5*Rx,ymin+0.9*Ry,title_str{mm},'interpreter','none','HorizontalAlignment','center','FontWeight','Bold','FontSize',10)
    end
    end
   saveas(hf,[dir_figure_out,figure_id,'_',sprintf('%03d',nfigure)],'png');
%   saveas(hf,[dir_figure_out,figure_id,'_',sprintf('%03d',nfigure)],'fig');

  end
end

mmax=k*Npage;
if (mmax<Nintf);
   hf= figure('units','normalized','outerposition',[0 0 1 1]);
   set(hf,'Visible','off')
    nfigure=nfigure+1;
    for kk=1:(Nintf-mmax);
    this_x=x{mmax+kk};
    this_y=y{mmax+kk};
    this_z=z{mmax+kk};
     subplot('position',psv{kk});
    hi=imagesc(this_x,this_y,this_z);
    colormap(jet);
    caxis([cmin cmax])
    set(gca,'XAxisLocation','top')
    set(hi,'alphadata',~isnan(this_z));
%    axis equal
%    axis tight


    if (kk==1);
     hcolor=colorbar('south','Position',[0.1+(Nx-1)*dx 0.915,dx-0.01,0.01]);
    end
    if (kk~=1);
     set(gca,'YTick',[])
     set(gca,'XTick',[])
    end
    set(gca,'YDir','normal')
    xmax=max(this_x);
    xmin=min(this_x);

   ymax=max(this_y);
   ymin=min(this_y);

   Rx=xmax-xmin;
   Ry=ymax-ymin;
   text(xmin+0.5*Rx,ymin+0.9*Ry,title_str{mmax+kk},'interpreter','none','HorizontalAlignment','center','FontWeight','Bold','FontSize',10)
 end
  saveas(hf,[dir_figure_out,figure_id,'_',sprintf('%03d',nfigure)],'png');
%  saveas(hf,[dir_figure_out,figure_id,'_',sprintf('%03d',nfigure)],'fig');
end

