function [xout,yout,zout,Npt,rms_out,xx1,xx2,yy1,yy2]=make_insar_downsample(xinsar,yinsar,zinsar,Nmin,Nres_min,Nres_max,method);
% downsample insar LOS measurements
%   
% Usage: [xout,yout,zout,Npt,rms_out,xx1,xx2,yy1,yy2]=make_insar_downsample(xinsar,yinsar,zinsar,Nmin,Nres_min,Nres_max,method);
%
%
%  Note: xinsar, yinsar are vectors, while zinsar is a matrix
%        xinsar=XINSASR(1,:);
%        yinsar=YINSAR(:,1);
%        [XINSAR,YINSAR]=meshgrid(xinsar,yinsar);

%  method: 'trend', or 'mean'
%  by Kang Wang on 08/27/2015
%  Last Updated by Kang Wang on 08/28/2015
%  Last Updated by Kang Wang on 09/01/2015
%  Last Updated by Kang Wang in Nov., 2015
%  Last Updated by Kang Wang on 10/04/2016

ismean=strcmp(method,'mean');
istrend=strcmp(method,'trend');

r1=20;
if (istrend);
 [xout,yout,zout,Npt,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_trend(xinsar,yinsar,zinsar,r1,Nres_min,Nres_max,[],[],[],[],[],[],[],[],[]);
elseif (ismean);
 [xout,yout,zout,Npt,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_mean(xinsar,yinsar,zinsar,r1,Nres_min,Nres_max,[],[],[],[],[],[],[],[],[]);
else
 error(['Not found method: ',method]);
end
Ndata=length(zout);

Nint=0;
while (Ndata<Nmin);
   N1=length(zout);
   r1=r1*0.85;
%   Ndata
   if (istrend);
    [xout,yout,zout,Npt,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_trend(xinsar,yinsar,zinsar,r1,Nres_min,Nres_max,[],[],[],[],[],[],[],[],[]);
   elseif(ismean);
    [xout,yout,zout,Npt,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_mean(xinsar,yinsar,zinsar,r1,Nres_min,Nres_max,[],[],[],[],[],[],[],[],[]);
   end
   Ndata=length(zout);
   N2=length(zout);  %after decreasing the threshold
   Nint=Nint+1;
   if (N2>0.8*Nmin&(N2-N1)<0.005*N1);
     break;     %stop the iteration if the the number of points are not increasing by 0.5 percent of the previous one
   end
end 

%Nint


function [xout,yout,zout,Ndata,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_trend(xin,yin,zin,threshold,Nres_min,Nres_max,...
xout_in,yout_in,zout_in,Ndata_in,rms_in,xx1_in,xx2_in,yy1_in,yy2_in);

%by Kang Wang on 08/27/2015
%Last Updated by Kang Wang on 08/28/2015

format long

%nx=length(unique(xin));
%ny=length(unique(yin));
nx=length(xin);
ny=length(yin);

xout=xout_in;
yout=yout_in;
zout=zout_in;
Ndata=Ndata_in;
rms_out=rms_in;
xx1=xx1_in;
xx2=xx2_in;
yy1=yy1_in;
yy2=yy2_in;

rms_default=10; %default rms
r_good_default=0.4; % threshold default percentage of good pixels
if (nx<=Nres_min|ny<=Nres_min);  %only continue if the size of the blocks are greater than Nres*Nres
  n_block=nx*ny;
  xout_block_total=mean(xin);
  yout_block_total=mean(yin);
  z_block_good=zin(~isnan(zin));
  N_block_good=length(z_block_good);
  r_good=N_block_good/n_block;
  if (N_block_good>0&r_good>r_good_default); 
    zout_block_total=mean(z_block_good);
    dz_block=z_block_good-zout_block_total;
    rms_block_total=sqrt(sum(sum(dz_block.^2))/N_block_good);
    if (rms_block_total<1.0e-6);
      rms_block_total=rms_default;
    end
    xout=[xout;xout_block_total];
    yout=[yout;yout_block_total];
    zout=[zout;zout_block_total];
    Ndata=[Ndata;N_block_good];
    rms_out=[rms_out;rms_block_total];
    xx1=[xx1;xin(1)];
    xx2=[xx2;xin(nx)];
    yy1=[yy1;yin(1)];
    yy2=[yy2;yin(ny)];
  end

else

  nx1=1;
  nx2=floor(nx/2)+1;
  nx3=nx2;
  nx4=nx;

  ny1=1;
  ny2=floor(ny/2)+1;
  ny3=ny2;
  ny4=ny;

  x1=xin(nx1:nx2);
  x2=xin(nx3:nx4);
  x3=x1;
  x4=x2;

  y1=yin(ny1:ny2);
  y2=y1;

  y3=yin(ny3:ny4);
  y4=y3;

  z1=zin(ny1:ny2,nx1:nx2);
  z2=zin(ny1:ny2,nx3:nx4);
  z3=zin(ny3:ny4,nx1:nx2);
  z4=zin(ny3:ny4,nx3:nx4);

  xmin1=xin(nx1);
  xmax1=xin(nx2);
  xmin2=xin(nx3);
  xmax2=xin(nx4);

  ymin1=yin(ny1);
  ymax1=yin(ny2);
  ymin2=yin(ny3);
  ymax2=yin(ny4);



  [rms1,N1,r1_good,x1_out,y1_out,z1_out]=rms_block_detrend(x1,y1,z1,Nres_min,Nres_max);
  [rms2,N2,r2_good,x2_out,y2_out,z2_out]=rms_block_detrend(x2,y2,z2,Nres_min,Nres_max);
  [rms3,N3,r3_good,x3_out,y3_out,z3_out]=rms_block_detrend(x3,y3,z3,Nres_min,Nres_max);
  [rms4,N4,r4_good,x4_out,y4_out,z4_out]=rms_block_detrend(x4,y4,z4,Nres_min,Nres_max);

  if (rms1<=threshold&N1>0&r1_good>r_good_default);
     xout=[xout;x1_out];
     yout=[yout;y1_out];
     zout=[zout;z1_out];

     zgood_this_block=z1(~isnan(z1));
     Ngood_this_block=length(zgood_this_block);
     dz_this_block=zgood_this_block-z1_out;
     rms_this_block=sqrt(sum(dz_this_block.^2)/Ngood_this_block);
     if (rms_this_block<1.0e-6);
       rms_this_block=rms_default;
     end
    
     rms_out=[rms_out;rms_this_block];
     Ndata=[Ndata;Ngood_this_block];
     xx1_this_block=xmin1;
     xx2_this_block=xmax1;
     yy1_this_block=ymin1;
     yy2_this_block=ymax1;

     xx1=[xx1;xx1_this_block];
     yy1=[yy1;yy1_this_block];
     xx2=[xx2;xx2_this_block];
     yy2=[yy2;yy2_this_block];
  elseif (rms1>threshold&N1>0);

     xout_in=xout;
     yout_in=yout;
     zout_in=zout;
     xx1_in=xx1;
     xx2_in=xx2;
     yy1_in=yy1;
     yy2_in=yy2;
     Ndata_in=Ndata;
     rms_in=rms_out;
     [xout,yout,zout,Ndata,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_trend(x1,y1,z1,threshold,Nres_min,Nres_max,...
     xout_in,yout_in,zout_in,Ndata_in,rms_in,xx1_in,xx2_in,yy1_in,yy2_in);

  end


  if (rms2<=threshold&N2>0&r2_good>r_good_default);
     xout=[xout;x2_out];
     yout=[yout;y2_out];
     zout=[zout;z2_out];

     zgood_this_block=z2(~isnan(z2));
     Ngood_this_block=length(zgood_this_block);
     dz_this_block=zgood_this_block-z2_out;
     rms_this_block=sqrt(sum(dz_this_block.^2)/Ngood_this_block);
     if (rms_this_block<1.0e-6);
        rms_this_block=rms_default;
     end

     rms_out=[rms_out;rms_this_block];
     Ndata=[Ndata;Ngood_this_block];
     xx1_this_block=xmin2;
     xx2_this_block=xmax2;
     yy1_this_block=ymin1;
     yy2_this_block=ymax1;

     xx1=[xx1;xx1_this_block];
     yy1=[yy1;yy1_this_block];
     xx2=[xx2;xx2_this_block];
     yy2=[yy2;yy2_this_block];
  elseif (rms2>threshold&N2>0);

     xout_in=xout;
     yout_in=yout;
     zout_in=zout;
     xx1_in=xx1;
     xx2_in=xx2;
     yy1_in=yy1;
     yy2_in=yy2;
     Ndata_in=Ndata;
     rms_in=rms_out;
     [xout,yout,zout,Ndata,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_trend(x2,y2,z2,threshold,Nres_min,Nres_max,...
     xout_in,yout_in,zout_in,Ndata_in,rms_in,xx1_in,xx2_in,yy1_in,yy2_in);
  end

  if (rms3<=threshold&N3>0&r3_good>r_good_default);
     xout=[xout;x3_out];
     yout=[yout;y3_out];
     zout=[zout;z3_out];

     zgood_this_block=z3(~isnan(z3));
     Ngood_this_block=length(zgood_this_block);
     dz_this_block=zgood_this_block-z3_out;
     rms_this_block=sqrt(sum(dz_this_block.^2)/Ngood_this_block);
     if (rms_this_block<1.0e-6);
        rms_this_block=rms_default;                                                         
     end
   
     rms_out=[rms_out;rms_this_block];                                                      
     Ndata=[Ndata;Ngood_this_block];                                                        
     xx1_this_block=xmin1;
     xx2_this_block=xmax1;
     yy1_this_block=ymin2;
     yy2_this_block=ymax2;
       
     xx1=[xx1;xx1_this_block];
     yy1=[yy1;yy1_this_block];
     xx2=[xx2;xx2_this_block];
     yy2=[yy2;yy2_this_block];
  elseif (rms3>threshold&N3>0);
      
     xout_in=xout;
     yout_in=yout;
     zout_in=zout;
     xx1_in=xx1;
     xx2_in=xx2;
     yy1_in=yy1;
     yy2_in=yy2;
     Ndata_in=Ndata;
     rms_in=rms_out;
     [xout,yout,zout,Ndata,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_trend(x3,y3,z3,threshold,Nres_min,Nres_max,...
     xout_in,yout_in,zout_in,Ndata_in,rms_in,xx1_in,xx2_in,yy1_in,yy2_in);
  end



  if (rms4<=threshold&N4>0&r4_good>r_good_default);
     xout=[xout;x4_out];
     yout=[yout;y4_out];
     zout=[zout;z4_out];

     zgood_this_block=z4(~isnan(z4));
     Ngood_this_block=length(zgood_this_block);
     dz_this_block=zgood_this_block-z4_out;
     rms_this_block=sqrt(sum(dz_this_block.^2)/Ngood_this_block);
     if (rms_this_block<1.0e-6);
       rms_this_block=rms_default;
     end

     rms_out=[rms_out;rms_this_block];
     Ndata=[Ndata;Ngood_this_block];
     xx1_this_block=xmin2;
     xx2_this_block=xmax2;
     yy1_this_block=ymin2;
     yy2_this_block=ymax2;

     xx1=[xx1;xx1_this_block];
     yy1=[yy1;yy1_this_block];
     xx2=[xx2;xx2_this_block];
     yy2=[yy2;yy2_this_block];
  elseif (rms4>threshold&N4>0);

     xout_in=xout;
     yout_in=yout;
     zout_in=zout;
     xx1_in=xx1;
     xx2_in=xx2;
     yy1_in=yy1;
     yy2_in=yy2;
     Ndata_in=Ndata;
     rms_in=rms_out;
    [xout,yout,zout,Ndata,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_trend(x4,y4,z4,threshold,Nres_min,Nres_max,...
     xout_in,yout_in,zout_in,Ndata_in,rms_in,xx1_in,xx2_in,yy1_in,yy2_in);
  end


end


function [rms_out,Ngood,r_good,xout,yout,zout]=rms_block_detrend(x,y,z,Nres_min,Nres_max);
format long
[xx,yy]=meshgrid(x,y);
indx_good=~isnan(z);
[nx,ny]=size(z);
n_block=nx*ny;

xdata=xx(indx_good);
ydata=yy(indx_good);
zdata=z(indx_good);

Ngood=length(zdata);
r_good=Ngood/n_block;

 if (Ngood>0);

  xout=mean(xdata);
  yout=mean(ydata);
  zout=mean(zdata);

  lx=length(unique(x));
  ly=length(unique(y)); 
%   lx=length(x);
%   ly=length(y);
  if (Ngood<=3|lx<=Nres_min|ly<=Nres_min);
     rms_out=0;
  elseif (Ngood>3&(lx>2&lx<=Nres_max)&(ly>2&ly<=Nres_max));
     xx=xdata;
     yy=ydata;
     zz=zdata;
     O=ones(Ngood,1);
     A=[xx,yy,O];
     C=pinv(A)*zz;
     zzfit=C(1)*xx+C(2)*yy+C(3);
     dz=zz-zzfit;
     rms_out=sqrt(sum(dz.^2)/Ngood);
%  elseif (length(unique(x))>30|length(unique(y))>30);
%     rms_out=1000;
  else
     rms_out=1000;
  end

 else 
   xout=NaN;    %no output
   yout=NaN;
   zout=NaN;
   rms_out=0;   
 end
%end   %end of the function



function [xout,yout,zout,Ndata,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_mean(xin,yin,zin,threshold,Nres_min,Nres_max,...
xout_in,yout_in,zout_in,Ndata_in,rms_in,xx1_in,xx2_in,yy1_in,yy2_in);

%by Kang Wang on 08/27/2015
%Last Updated by Kang Wang on 08/28/2015

format long

%nx=length(unique(xin));
%ny=length(unique(yin));
nx=length(xin);
ny=length(yin);

xout=xout_in;
yout=yout_in;
zout=zout_in;
Ndata=Ndata_in;
rms_out=rms_in;
xx1=xx1_in;
xx2=xx2_in;
yy1=yy1_in;
yy2=yy2_in;

rms_default=10; %default rms
r_good_default=0.4; %defualt percentage of good pixels
if (nx<=Nres_min|ny<=Nres_min);  %only continue if the size of the blocks are greater than Nres*Nres
  xout_block_total=mean(xin);
  yout_block_total=mean(yin);
  z_block_good=zin(~isnan(zin));
  N_block_good=length(z_block_good);
  n_block=nx*ny;
  r_good=N_block_good/n_block;
  if (N_block_good>0&r_good>r_good_default);
    zout_block_total=mean(z_block_good);
    dz_block=z_block_good-zout_block_total;
    rms_block_total=sqrt(sum(sum(dz_block.^2))/N_block_good);
    if (rms_block_total<1.0e-6);
      rms_block_total=10;
    end
    xout=[xout;xout_block_total];
    yout=[yout;yout_block_total];
    zout=[zout;zout_block_total];
    Ndata=[Ndata;N_block_good];
    rms_out=[rms_out;rms_block_total];
    xx1=[xx1;xin(1)];
    xx2=[xx2;xin(nx)];
    yy1=[yy1;yin(1)];
    yy2=[yy2;yin(ny)];
  end

else

  nx1=1;
  nx2=floor(nx/2)+1;
  nx3=nx2;
  nx4=nx;

  ny1=1;
  ny2=floor(ny/2)+1;
  ny3=ny2;
  ny4=ny;

  x1=xin(nx1:nx2);
  x2=xin(nx3:nx4);
  x3=x1;
  x4=x2;

  y1=yin(ny1:ny2);
  y2=y1;

  y3=yin(ny3:ny4);
  y4=y3;

  z1=zin(ny1:ny2,nx1:nx2);
  z2=zin(ny1:ny2,nx3:nx4);
  z3=zin(ny3:ny4,nx1:nx2);
  z4=zin(ny3:ny4,nx3:nx4);

  xmin1=xin(nx1);
  xmax1=xin(nx2);
  xmin2=xin(nx3);
  xmax2=xin(nx4);

  ymin1=yin(ny1);
  ymax1=yin(ny2);
  ymin2=yin(ny3);
  ymax2=yin(ny4);



  [rms1,N1,r1_good,x1_out,y1_out,z1_out]=rms_block_demean(x1,y1,z1,Nres_min,Nres_max);
  [rms2,N2,r2_good,x2_out,y2_out,z2_out]=rms_block_demean(x2,y2,z2,Nres_min,Nres_max);
  [rms3,N3,r3_good,x3_out,y3_out,z3_out]=rms_block_demean(x3,y3,z3,Nres_min,Nres_max);
  [rms4,N4,r4_good,x4_out,y4_out,z4_out]=rms_block_demean(x4,y4,z4,Nres_min,Nres_max);

  if (rms1<=threshold&N1>0&r1_good>r_good_default);
     xout=[xout;x1_out];
     yout=[yout;y1_out];
     zout=[zout;z1_out];

     zgood_this_block=z1(~isnan(z1));
     Ngood_this_block=length(zgood_this_block);
     dz_this_block=zgood_this_block-z1_out;
     rms_this_block=sqrt(sum(dz_this_block.^2)/Ngood_this_block);
     if (rms_this_block<1.0e-6);
       rms_this_block=rms_default;
     end
    
     rms_out=[rms_out;rms_this_block];
     Ndata=[Ndata;Ngood_this_block];
     xx1_this_block=xmin1;
     xx2_this_block=xmax1;
     yy1_this_block=ymin1;
     yy2_this_block=ymax1;

     xx1=[xx1;xx1_this_block];
     yy1=[yy1;yy1_this_block];
     xx2=[xx2;xx2_this_block];
     yy2=[yy2;yy2_this_block];
  elseif (rms1>threshold&N1>0);

     xout_in=xout;
     yout_in=yout;
     zout_in=zout;
     xx1_in=xx1;
     xx2_in=xx2;
     yy1_in=yy1;
     yy2_in=yy2;
     Ndata_in=Ndata;
     rms_in=rms_out;
     [xout,yout,zout,Ndata,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_mean(x1,y1,z1,threshold,Nres_min,Nres_max,...
     xout_in,yout_in,zout_in,Ndata_in,rms_in,xx1_in,xx2_in,yy1_in,yy2_in);

  end


  if (rms2<=threshold&N2>0&r2_good>r_good_default);
     xout=[xout;x2_out];
     yout=[yout;y2_out];
     zout=[zout;z2_out];

     zgood_this_block=z2(~isnan(z2));
     Ngood_this_block=length(zgood_this_block);
     dz_this_block=zgood_this_block-z2_out;
     rms_this_block=sqrt(sum(dz_this_block.^2)/Ngood_this_block);
     if (rms_this_block<1.0e-6);
        rms_this_block=rms_default;
     end

     rms_out=[rms_out;rms_this_block];
     Ndata=[Ndata;Ngood_this_block];
     xx1_this_block=xmin2;
     xx2_this_block=xmax2;
     yy1_this_block=ymin1;
     yy2_this_block=ymax1;

     xx1=[xx1;xx1_this_block];
     yy1=[yy1;yy1_this_block];
     xx2=[xx2;xx2_this_block];
     yy2=[yy2;yy2_this_block];
  elseif (rms2>threshold&N2>0);

     xout_in=xout;
     yout_in=yout;
     zout_in=zout;
     xx1_in=xx1;
     xx2_in=xx2;
     yy1_in=yy1;
     yy2_in=yy2;
     Ndata_in=Ndata;
     rms_in=rms_out;
     [xout,yout,zout,Ndata,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_mean(x2,y2,z2,threshold,Nres_min,Nres_max,...
     xout_in,yout_in,zout_in,Ndata_in,rms_in,xx1_in,xx2_in,yy1_in,yy2_in);
  end

  if (rms3<=threshold&N3>0&r3_good>r_good_default);
     xout=[xout;x3_out];
     yout=[yout;y3_out];
     zout=[zout;z3_out];

     zgood_this_block=z3(~isnan(z3));
     Ngood_this_block=length(zgood_this_block);
     dz_this_block=zgood_this_block-z3_out;
     rms_this_block=sqrt(sum(dz_this_block.^2)/Ngood_this_block);
     if (rms_this_block<1.0e-6);
        rms_this_block=rms_default;                                                         
     end
   
     rms_out=[rms_out;rms_this_block];                                                      
     Ndata=[Ndata;Ngood_this_block];                                                        
     xx1_this_block=xmin1;
     xx2_this_block=xmax1;
     yy1_this_block=ymin2;
     yy2_this_block=ymax2;
       
     xx1=[xx1;xx1_this_block];
     yy1=[yy1;yy1_this_block];
     xx2=[xx2;xx2_this_block];
     yy2=[yy2;yy2_this_block];
  elseif (rms3>threshold&N3>0);
      
     xout_in=xout;
     yout_in=yout;
     zout_in=zout;
     xx1_in=xx1;
     xx2_in=xx2;
     yy1_in=yy1;
     yy2_in=yy2;
     Ndata_in=Ndata;
     rms_in=rms_out;
     [xout,yout,zout,Ndata,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_mean(x3,y3,z3,threshold,Nres_min,Nres_max,...
     xout_in,yout_in,zout_in,Ndata_in,rms_in,xx1_in,xx2_in,yy1_in,yy2_in);
  end



  if (rms4<=threshold&N4>0&r4_good>r_good_default);
     xout=[xout;x4_out];
     yout=[yout;y4_out];
     zout=[zout;z4_out];

     zgood_this_block=z4(~isnan(z4));
     Ngood_this_block=length(zgood_this_block);
     dz_this_block=zgood_this_block-z4_out;
     rms_this_block=sqrt(sum(dz_this_block.^2)/Ngood_this_block);
     if (rms_this_block<1.0e-6);
       rms_this_block=rms_default;
     end

     rms_out=[rms_out;rms_this_block];
     Ndata=[Ndata;Ngood_this_block];
     xx1_this_block=xmin2;
     xx2_this_block=xmax2;
     yy1_this_block=ymin2;
     yy2_this_block=ymax2;

     xx1=[xx1;xx1_this_block];
     yy1=[yy1;yy1_this_block];
     xx2=[xx2;xx2_this_block];
     yy2=[yy2;yy2_this_block];
  elseif (rms4>threshold&N4>0);

     xout_in=xout;
     yout_in=yout;
     zout_in=zout;
     xx1_in=xx1;
     xx2_in=xx2;
     yy1_in=yy1;
     yy2_in=yy2;
     Ndata_in=Ndata;
     rms_in=rms_out;
    [xout,yout,zout,Ndata,rms_out,xx1,xx2,yy1,yy2]=quad_decomp_mean(x4,y4,z4,threshold,Nres_min,Nres_max,...
     xout_in,yout_in,zout_in,Ndata_in,rms_in,xx1_in,xx2_in,yy1_in,yy2_in);
  end
end


function [rms_out,Ngood,r_good,xout,yout,zout]=rms_block_demean(x,y,z,Nres_min,Nres_max);
format long
[xx,yy]=meshgrid(x,y);
indx_good=~isnan(z);
[nx,ny]=size(z);
n_block=nx*ny;
xdata=xx(indx_good);
ydata=yy(indx_good);
zdata=z(indx_good);

Ngood=length(zdata);
r_good=Ngood/n_block;
 if (Ngood>0);

  xout=mean(xdata);
  yout=mean(ydata);
  zout=mean(zdata);
  lx=length(unique(x));
  ly=length(unique(y));
%   lx=length(x);
%   ly=length(y);
  if (Ngood<=3|lx<=Nres_min|ly<=Nres_min);
    rms_out=0;
  elseif (Ngood>5&(lx>2&lx<Nres_max)&(ly>2&ly<Nres_max));
     zz=zdata;
     zzfit=mean(zz); 
     dz=zz-zzfit;
     rms_out=sqrt(sum(dz.^2)/Ngood);
  else
     rms_out=1000;
  end

 else 
   xout=NaN;    %no output
   yout=NaN;
   zout=NaN;
   rms_out=0;   
 end
%end   %end of the function



