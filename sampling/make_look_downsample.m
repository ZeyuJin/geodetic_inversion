function [xout,yout,zout]=make_look_downsample(xlook,ylook,zlook,xin,yin,xx1,xx2,yy1,yy2);
% downsample the looking vectors or other observables using the given boundings
%
%  Usage: [xout,yout,zout]=make_look_downsample(xlook,ylook,zlook,xin,yin,xx1,xx2,yy1,yy2);
%         xout=xin;
%         yout=yin;
%  
%  by Kang Wang on 08/27/2015
 
format long
Nblocks=length(xin); 

xout=[];
yout=[];
zout=[];
for k=1:Nblocks; 
  indx_x=find(xlook>=xx1(k)&xlook<=xx2(k));
  indx_y=find(ylook>=yy1(k)&ylook<=yy2(k));
  zfind=zlook(indx_y,indx_x);
  indx_good=~isnan(zfind);
  zgood=zfind(indx_good);
  Ngood=length(zgood);
  if (Ngood>0);
    z_this_block=mean(zgood);
  else 
    z_this_block=NaN;
  end
    zout=[zout;z_this_block];
    xout=[xout;xin(k)];
    yout=[yout;yin(k)];
end


