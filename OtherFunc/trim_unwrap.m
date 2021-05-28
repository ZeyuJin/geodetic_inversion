function [x,y,znew]=trim_unwrap(x,y,z,varargin);
%  
%  [x,y,znew]=trim_unwrap(x,y,z);
%
%
[X,Y]=meshgrid(x,y);
indx_good=~isnan(z);
xgood=X(indx_good);
ygood=Y(indx_good);
xmin=min(xgood);
xmax=max(xgood);
Rx=xmax-xmin;
ymin=min(ygood);
ymax=max(ygood);
Ry=ymax-ymin;

xxmin=xmin-0.03*Rx;
xxmax=xmax+0.03*Rx;

yymin=ymin-0.03*Ry;
yymax=ymax+0.03*Ry;

if (length(varargin)==0)
  region = [xxmin xxmax yymin yymax];
else
  region =[varargin{1}];
end

h1=figure('units','normalized','outerposition',[0 0 1 1]);
psv1=[0.05 0 0.45 1];
subplot('position',psv1);
h=imagesc(x,y,z);
set(gca,'YDir','Normal')
axis equal
colorbar
set(h,'alphadata',~isnan(z));
colormap('jet');
axis(region);

znew=z;
flag='y';
while (strcmp(flag,'y'));
 [xpoly,ypoly]=get_polygon(h1,6);
 indx_mask=inpolygon(X,Y,xpoly,ypoly);
 znew(indx_mask)=NaN;
 psv2=[0.55 0 0.45 1];
 subplot('position',psv2);
 h=imagesc(x,y,znew);
 set(gca,'YDir','Normal')
 axis equal
 colorbar
 set(h,'alphadata',~isnan(znew));
 colormap('jet');
 axis(region);

 prompt = 'Continue (y/n)?\n';
 flag = input(prompt,'s');
end

%grdwrite2(x,y,znew,['trim_',grd_input]);
