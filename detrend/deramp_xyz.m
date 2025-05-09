 function [ramp,cffs]=deramp_xyz(Ydata,xd,yd,shift);
%function [Ydata,cffs]=deramp_xyz(Ydata,xd,yd,shift);

max_amp=200; % threshould for "noisy" pixels
[ymax,xmax] = size(Ydata);
if (length(xd(:)) == 0)
 [xd,yd]=meshgrid(1:xmax,1:ymax);
end
cffs=[];
%goodR=find(isnan(Ydata)==0);
goodR=~isnan(Ydata);
data=Ydata(goodR); 
x=xd(goodR); 
y=yd(goodR);
%eval(['cd ' dpath]);
%load([dpath fname]);
switch shift % de-ramp stuff
 case 89      % no change
  return
 case 91      % constant shift
  A         = [ones(size(x))];
  cffs      = A\data;
%  fprintf('DC shift:\n'); 
%  fprintf('%e \n',cffs(1));
  ramp=cffs(1);
  Ydata=Ydata-cffs(1); 
 case 93      % planar ramp
  A         = [x, y, ones(size(x))];
  cffs      = A\data;
%  fprintf('plane coefficients:\n'); 
%  fprintf('%e %e %e \n',cffs(1),cffs(2),cffs(3));
  Ydata=Ydata-(cffs(1)*xd+cffs(2)*yd+cffs(3));
  ramp=(cffs(1)*xd+cffs(2)*yd+cffs(3));
 case 94      % bilinear ramp
  A         = [x.*y, x, y, ones(size(x))];
  cffs      = A\data;
%  fprintf('bilinear coefficients:\n'); 
%  fprintf('%e %e %e %e \n',cffs(1),cffs(2),cffs(3),cffs(4));
  Ydata=Ydata-(cffs(1)*xd.*yd+cffs(2)*xd+cffs(3)*yd+cffs(4));
  ramp=(cffs(1)*xd.*yd+cffs(2)*xd+cffs(3)*yd+cffs(4));
 case 95      % quadratic ramp
  A         = [x.^2 y.^2 x.*y x, y, ones(size(x))];
  cffs      = A\data;
  fprintf('quad coefficients:\n'); 
  fprintf('%e %e %e %e %e %e \n',cffs(1),cffs(2),cffs(3),cffs(4),cffs(5),cffs(6));
  Ydata=Ydata-(cffs(1)*xd.^2+cffs(2)*yd.^2+cffs(3)*xd.*yd+cffs(4)*xd+cffs(5)*yd+cffs(6));
  ramp=(cffs(1)*xd.^2+cffs(2)*yd.^2+cffs(3)*xd.*yd+cffs(4)*xd+cffs(5)*yd+cffs(6));
% otherwise
%  data=Ydata-shift;        % add DC shift
end   %switch

return 

% detrend second time, tossing out high amplitude ("noisy"?) pixels
cffs=[];
Ddata=Ydata;
bad=find(abs(Ddata)>max_amp);
Ddata(bad)=nan;
goodR=find(isnan(Ddata)==0);
data=Ddata(goodR); x=xd(goodR); y=yd(goodR);
%eval(['cd ' dpath]);
%load([dpath fname]);
switch shift % de-ramp stuff
 case 89      % no change
  return
 case 91      % constant shift
  A         = [ones(size(x))];
  cffs      = A\data;
%  fprintf('DC shift:\n'); 
%  fprintf('%e \n',cffs(1));
  Ydata=Ydata-cffs(1); 
  ramp=cffs(1);
 case 93      % planar ramp
  A         = [x, y, ones(size(x))];
  cffs      = A\data;
%  fprintf('plane coefficients:\n'); 
%  fprintf('%e %e %e \n',cffs(1),cffs(2),cffs(3));
  Ydata=Ydata-(cffs(1)*xd+cffs(2)*yd+cffs(3));
  ramp=Ydata-(cffs(1)*xd+cffs(2)*yd+cffs(3));
 case 94      % bilinear ramp
  A         = [x.*y, x, y, ones(size(x))];
  cffs      = A\data;
%  fprintf('bilinear coefficients:\n'); 
%  fprintf('%e %e %e %e \n',cffs(1),cffs(2),cffs(3),cffs(4));
  Ydata=Ydata-(cffs(1)*xd.*yd+cffs(2)*xd+cffs(3)*yd+cffs(4));
  ramp=(cffs(1)*xd.*yd+cffs(2)*xd+cffs(3)*yd+cffs(4));
 case 95      % quadratic ramp
  A         = [x.^2 y.^2 x.*y x, y, ones(size(x))];
  cffs      = A\data;
  fprintf('quad coefficients:\n'); 
  fprintf('%e %e %e %e %e %e \n',cffs(1),cffs(2),cffs(3),cffs(4),cffs(5),cffs(6));

  Ydata=Ydata-(cffs(1)*xd.^2+cffs(2)*yd.^2+cffs(3)*xd.*yd+cffs(4)*xd+cffs(5)*yd+cffs(6));
  ramp=(cffs(1)*xd.^2+cffs(2)*yd.^2+cffs(3)*xd.*yd+cffs(4)*xd+cffs(5)*yd+cffs(6));
% otherwise
%  data=Ydata-shift;        % add DC shift
end   %switch


%fid = fopen(fname,'w');
%status = fclose(fid);
%save(fname,'cffs');
%f=1;

