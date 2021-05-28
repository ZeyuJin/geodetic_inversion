function plot_insar_grd(x,y,varargin);
% plot out the insar observations 
% Usage: plot_insar_grd(x,y,varargin);
% x,y ----- vectors of coordinates
% z   ----- 2D matrix of the observation/model/res
%
% example: plot_insar_grd(x,y,zobs,zmodel,zres);
%
% by Kang Wang on 10/28/2015
%
set(0,'defaultAxesFontName', 'Helvetica')
set(0,'defaultAxesFontSize', 15)
nfields = length(varargin);
if (nfields==1);
  z1=varargin{1};
  zmin=min(min(z1(~isnan(z1))));
  zmax=max(max(z1(~isnan(z1))));
  rz=zmax-zmin;

  cmin=zmin+0.05*rz;
  cmax=zmax-0.05*rz;
  
  hf=figure('units','normalized','outerposition',[0 0 1 1]);
  psv1=[0.1 0.2 0.7 0.7];
  subplot('position',psv1);
  h1=imagesc(x,y,z1);
  set(gca,'YDir','Normal')
  axis equal
  colormap(jet)
  set(h1,'alphadata',~isnan(z1));
  hcolor=colorbar('south','position',[0.2 0.15 0.3 0.01]);
  caxis([cmin cmax])
end


if (nfields==3);
    z1=varargin{1};
    z2=varargin{2};
    z3=varargin{3};
    
    zmin=min(min(z1(~isnan(z1))));
    zmax=max(max(z1(~isnan(z1))));
    rz=zmax-zmin;

    cmin=zmin+0.05*rz;
    cmax=zmax-0.05*rz;

    zmin2=min(min(z3(~isnan(z3))));
    zmax2=max(max(z3(~isnan(z3))));
    rz2=zmax2-zmin2;

    cmin2=zmin2+0.05*rz2;
    cmax2=zmax2-0.05*rz2;    

    hf=figure('units','normalized','outerposition',[0 0 1 1]);
    psv1=[0.02 0.2 0.3 0.6];
    subplot('position',psv1);
    h1=imagesc(x,y,z1);
    set(gca,'YDir','Normal')
    axis equal
    %axis tight
    colormap(jet)
    set(h1,'alphadata',~isnan(z1))
    caxis([cmin cmax])

    psv2=[0.33 0.2 0.3 0.6];
    subplot('position',psv2)
    h2=imagesc(x,y,z2);
    set(gca,'YDir','Normal');
    axis equal
    %axis tight
    colormap(jet)
    set(gca,'Ytick',[])
    set(h2,'alphadata',~isnan(z2))
    caxis([cmin cmax])    

    hcolor=colorbar('south','position',[0.15 0.15 0.3 0.01]);
    
    psv3=[0.65 0.2 0.3 0.6];
    subplot('position',psv3);
    h3=imagesc(x,y,z3);
    set(gca,'YDir','Normal')
    axis equal
    colormap(jet)
    set(gca,'Ytick',[])
    set(h3,'alphadata',~isnan(z3))
    
    hcolor=colorbar('south','position',[0.7 0.15 0.2 0.01]);
    caxis([cmin2 cmax2])

end
