function [rnew,znew,sigma,L]=get_insar_varigram(xx,yy,zz,rmin,rmax,dr,Nmax);
%  get the insar semi-varigram assuming the noise structure is isotropic,
%
%  and fit it as exponental function
%
%  Usage: [rout,zout,sigma,L]=get_insar_varigram(xin,yin,zin,rmin,rmax,dr,Nmax);
%
%     input:   xin ---- x coordinates of the noise [m]
%              yin ---- y coordinates of the noise [m]
%              zin ---- noise [cm]
%              rmin ---- minimum distance
%              rmax ---- maximum distance
%              dr   ---- distance step
%              Nmax ---- number pairs to use
% 
%    C(r) = sigma * exp(-r/L);
%
%    by Kang Wang in  Feb. 2018

indx_good=~isnan(zz);
xin=xx(indx_good);
yin=yy(indx_good);
zin=zz(indx_good);

var_all=var(zin);

npt=length(zin);
indx=floor(npt*rand(Nmax,2))+1;
indx_new=unique(indx,'rows');
indx1=indx_new(:,1);
indx2=indx_new(:,2);

z1=zin(indx1);
z2=zin(indx2);
zdiff = zin(indx1)-zin(indx2);

x1=xin(indx1);
y1=yin(indx1);

x2=xin(indx2);
y2=yin(indx2);

r=sqrt((x2-x1).^2+(y2-y1).^2);

rin=[rmin:dr:rmax];
nr=length(rin);

rout=NaN(nr-1,1);
zout=NaN(nr-1,1);
for k=1:nr-1;
  r1 = rin(k);
  r2 = rin(k+1);
  rout(k)=mean([r1,r2]); % take the middle point
  indx_bin=find(r>=r1 & r<r2);
  N=length(indx_bin);
  if (N>2);  % do this for bins with at least 3 pair;
     zout(k)=var_all-var(zdiff(indx_bin))/2;
  end
  
end


zout=zout(~isnan(zout));
rout=rout(~isnan(zout));

[rnew,indx_sort]=sort(rout);
znew=zout(indx_sort);

F= @(x, xdata) x(1)*exp(-xdata/x(2));
x0=[1,10e3];
[x,resnorm,~,exitflag,output]=lsqcurvefit(F,x0,rnew,znew);
sigma=x(1);
L=x(2);
