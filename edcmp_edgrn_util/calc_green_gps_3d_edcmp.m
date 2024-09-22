function G=calc_green_gps_3d_edcmp(slip_model,data_insar,finp);
% Calculate the Green's function for InSAR observations using EDGRN
%
% Usage: G=calc_green_insar(slip_model,data_insar,fin);
%
% by Kang Wang on 08/21/2019
% Last Updated by Kang Wang on 02/20/2020

format long
d2r=pi/180;

% if not exist the MAT file
% then run EDGRN first
if ~isfile(['edgrn_ss_',finp,'.mat'])
    getedgrn(finp);
end

% fix this "length" bug
npatch=size(slip_model, 1); %number of patches
npara=2*npatch; %number of parameters (strike-slip + dip-slip)

% for InSAR data
yrec=data_insar(:,1);
xrec=data_insar(:,2);
ve=data_insar(:,4);
vn=data_insar(:,5);
vz=data_insar(:,6);

nobs = length(yrec);
ZSEPS=1.0d-2;

% for GPS data
% nstn=length(xe_gps);
% nobs=3*nstn;

ss=load(['edgrn_ss_',finp,'.mat']);
ds=load(['edgrn_ds_',finp,'.mat']);
cl=load(['edgrn_cl_',finp,'.mat']);
nr=ss.nr;
r1=ss.r1;
r2=ss.r2;
nz=ss.nz;
z1=ss.z1;
z2=ss.z2;
zrec0=ss.zrec0;
lambda=ss.lambda;
mu=ss.mu;


% strike-slip component
ssu1=zeros(nr+2,nz+2);
ssu2=zeros(nr+2,nz+2);
ssu3=zeros(nr+2,nz+2);
ssu1(1:nr,1:nz)=reshape(ss.uz,nr,nz);
ssu2(1:nr,1:nz)=reshape(ss.ur,nr,nz);
ssu3(1:nr,1:nz)=reshape(ss.ut,nr,nz);

% dip-slip component
dsu1=zeros(nr+2,nz+2);
dsu2=zeros(nr+2,nz+2);
dsu3=zeros(nr+2,nz+2);
dsu1(1:nr,1:nz)=reshape(ds.uz,nr,nz);
dsu2(1:nr,1:nz)=reshape(ds.ur,nr,nz);
dsu3(1:nr,1:nz)=reshape(ds.ut,nr,nz);

% clvd
clu1=zeros(nr+2,nz+2);
clu2=zeros(nr+2,nz+2);
clu3=zeros(nr+2,nz+2);
clu1(1:nr,1:nz)=reshape(cl.uz,nr,nz);
clu2(1:nr,1:nz)=reshape(cl.ur,nr,nz);
clu3(1:nr,1:nz)=reshape(cl.ut,nr,nz);


dr=(r2-r1)/(nr-1);
dz=(z2-z1)/(nz-1);

if ((zrec0>=z1) &(zrec0<=z2) &(r1==0))
    iz=floor(zrec0-z1)/dz+1;
    dzs=(zrec0-(z1+dz*(iz-1)))/dz;
    if (abs(dzs)<=ZSEPS)
      ssu1(1,iz)=ssu1(2,iz);
      ssu2(1,iz)=ssu2(2,iz);
      ssu3(1,iz)=ssu3(2,iz);
      
      dsu1(1,iz)=dsu1(2,iz);
      dsu2(1,iz)=dsu2(2,iz);
      dsu3(1,iz)=dsu3(2,iz);
      
      clu1(1,iz)=clu1(2,iz);
      clu2(1,iz)=clu2(2,iz);

    end
end

G=zeros(nobs,npara);

% fprintf('%s\n',['working on the strike-slip component']);

for ipatch = 1:npatch;
   if mod(ipatch, 400) == 0
       disp(['... Green-GPS strike-slip patch: ',' [',num2str(ipatch),'/',num2str(npatch),']']);
   end
   xs_patch=slip_model(ipatch,5);
   ys_patch=slip_model(ipatch,4);
   zs_patch=-slip_model(ipatch,6);  % reverse the z-axis
   lp_patch=slip_model(ipatch,7);
   wp_patch=slip_model(ipatch,8);
   strike_patch=slip_model(ipatch,9);
   dip_patch=slip_model(ipatch,10);

   slip_patch=1;
   rake_patch=0;
   [pxs,pys,pzs,pmoment]=edcdisc(xs_patch,ys_patch,zs_patch,slip_patch,lp_patch,...
   wp_patch,strike_patch,dip_patch,rake_patch,nz,z1,z2,dr,dz);
   nrec=length(xrec); %number of observation points
   nps=length(pxs);   %number of point sources

   %superposition of all discrete point sources
   Ue=zeros(nrec,1);
   Un=zeros(nrec,1);
   Up=zeros(nrec,1);

   for k=1:nrec;
   
      dis=sqrt((xrec(k)-pxs).^2+(yrec(k)-pys).^2);
    if (dis<0);
      azi=0;
    else
      azi=atan2(yrec(k)-pys,xrec(k)-pxs);  
    end
   
     iz=floor((pzs-z1)/dz)+1;
     dzs=(pzs-(z1+dz*(iz-1)))/dz;
   
     if (dis<r1);
      idis=1;
      ddis=0;
     else
       idis=floor((dis-r1)/dr)+1;
       ddis=(dis-(r1+dr*(idis-1)))/dr;
     end
   
   
      w00=(1.d0-ddis).*(1.d0-dzs);
      w10=ddis.*(1.d0-dzs);
      w01=(1.d0-ddis).*dzs;
      w11=ddis.*dzs;
      co=cos(azi);
      si=sin(azi);
      co2=cos(2d0*azi);
      si2=sin(2d0*azi);
      
%    contributions from the strike-slip components
      ps=pmoment(1,:)'.*si2+pmoment(4,:)'.*co2;
      sh=pmoment(1,:)'.*co2-pmoment(4,:)'.*si2;
     
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);

      for nn=1:nps;
        AA1(nn)=ssu1(idis(nn),iz(nn));
        AA2(nn)=ssu1(idis(nn)+1,iz(nn));
        AA3(nn)=ssu1(idis(nn),iz(nn)+1);
        AA4(nn)=ssu1(idis(nn)+1,iz(nn)+1);
      end
      
      
      suz=ps.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);

      for nn=1:nps;
        AA1(nn)=ssu2(idis(nn),iz(nn));
        AA2(nn)=ssu2(idis(nn)+1,iz(nn));
        AA3(nn)=ssu2(idis(nn),iz(nn)+1);
        AA4(nn)=ssu2(idis(nn)+1,iz(nn)+1);
      end
      sur=ps.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);

      for nn=1:nps;
        AA1(nn)=ssu3(idis(nn),iz(nn));
        AA2(nn)=ssu3(idis(nn)+1,iz(nn));
        AA3(nn)=ssu3(idis(nn),iz(nn)+1);
        AA4(nn)=ssu3(idis(nn)+1,iz(nn)+1);
      end
      sut=sh.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);

%    
%   contribution from the dip-slip comoents
      ps=pmoment(2,:)'.*co+pmoment(5,:)'.*si;
      sh=pmoment(2,:)'.*si-pmoment(5,:)'.*co;
      
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);

      for nn=1:nps;
        AA1(nn)=dsu1(idis(nn),iz(nn));
        AA2(nn)=dsu1(idis(nn)+1,iz(nn));
        AA3(nn)=dsu1(idis(nn),iz(nn)+1);
        AA4(nn)=dsu1(idis(nn)+1,iz(nn)+1);
      end
      duz=ps.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);

      for nn=1:nps;
        AA1(nn)=dsu2(idis(nn),iz(nn));
        AA2(nn)=dsu2(idis(nn)+1,iz(nn));
        AA3(nn)=dsu2(idis(nn),iz(nn)+1);
        AA4(nn)=dsu2(idis(nn)+1,iz(nn)+1);
      end
      dur=ps.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);
      for nn=1:nps;
        AA1(nn)=dsu3(idis(nn),iz(nn));
        AA2(nn)=dsu3(idis(nn)+1,iz(nn));
        AA3(nn)=dsu3(idis(nn),iz(nn)+1);
        AA4(nn)=dsu3(idis(nn)+1,iz(nn)+1);
      end
      dut=sh.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
%  contributions from the clvd component
      ps=pmoment(3,:)';
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);

      for nn=1:nps;
        AA1(nn)=clu1(idis(nn),iz(nn));
        AA2(nn)=clu1(idis(nn)+1,iz(nn));
        AA3(nn)=clu1(idis(nn),iz(nn)+1);
        AA4(nn)=clu1(idis(nn)+1,iz(nn)+1);
      end
      cluz=ps.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);
      for nn=1:nps;
        AA1(nn)=clu2(idis(nn),iz(nn));
        AA2(nn)=clu2(idis(nn)+1,iz(nn));
        AA3(nn)=clu2(idis(nn),iz(nn)+1);
        AA4(nn)=clu2(idis(nn)+1,iz(nn)+1);
      end
      clur=ps.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
      Uz=suz+duz+cluz;
      Ur=sur+dur+clur;
      Ut=sut+dut;
    
      Ux=sum(Ur.*co-Ut.*si);
      Uy=sum(Ur.*si+Ut.*co);
      Uz=sum(Uz);
      
      Un(k)=Ux;
      Ue(k)=Uy;
      Up(k)=-Uz;
   end
  
   % maybe there is one bug here
%    G(:,ipatch+npatch)=[Ue;Un;Up];
   G(:,ipatch) = Ue.*ve + Un.*vn + Up.*vz;

end

% fprintf('%s\n',['working on the dip-slip component']);

for ipatch = 1:npatch;
   if mod(ipatch+1, 400) == 0
       disp(['... Green-GPS dip-slip patch: ',' [',num2str(ipatch),'/',num2str(npatch),']']);
   end

   xs_patch=slip_model(ipatch,5);
   ys_patch=slip_model(ipatch,4);
   zs_patch=-slip_model(ipatch,6);
   lp_patch=slip_model(ipatch,7);
   wp_patch=slip_model(ipatch,8);
   strike_patch=slip_model(ipatch,9);
   dip_patch=slip_model(ipatch,10);
   slip_patch=1;
   rake_patch=90;
   [pxs,pys,pzs,pmoment]=edcdisc(xs_patch,ys_patch,zs_patch,slip_patch,lp_patch,...
    wp_patch,strike_patch,dip_patch,rake_patch,nz,z1,z2,dr,dz);
    nrec=length(xrec); %number of observation points
    nps=length(pxs);   %number of point sources

  %superposition of all discrete point sources

    Ue=zeros(nrec,1);
    Un=zeros(nrec,1);
    Up=zeros(nrec,1);

   for k=1:nrec;
   
      dis=sqrt((xrec(k)-pxs).^2+(yrec(k)-pys).^2);
    if (dis<0);
      azi=0;
    else
      azi=atan2(yrec(k)-pys,xrec(k)-pxs);  
    end
   
     iz=floor((pzs-z1)/dz)+1;
     dzs=(pzs-(z1+dz*(iz-1)))/dz;
   
     if (dis<r1);
      idis=1;
      ddis=0;
     else
       idis=floor((dis-r1)/dr)+1;
       ddis=(dis-(r1+dr*(idis-1)))/dr;
     end
   
   
      w00=(1.d0-ddis).*(1.d0-dzs);
      w10=ddis.*(1.d0-dzs);
      w01=(1.d0-ddis).*dzs;
      w11=ddis.*dzs;
      co=cos(azi);
      si=sin(azi);
      co2=cos(2d0*azi);
      si2=sin(2d0*azi);
      
%    contributions from the strike-slip components
      ps=pmoment(1,:)'.*si2+pmoment(4,:)'.*co2;
      sh=pmoment(1,:)'.*co2-pmoment(4,:)'.*si2;
     
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);

      for nn=1:nps;
        AA1(nn)=ssu1(idis(nn),iz(nn));
        AA2(nn)=ssu1(idis(nn)+1,iz(nn));
        AA3(nn)=ssu1(idis(nn),iz(nn)+1);
        AA4(nn)=ssu1(idis(nn)+1,iz(nn)+1);
      end
      
      
      suz=ps.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);

      for nn=1:nps;
        AA1(nn)=ssu2(idis(nn),iz(nn));
        AA2(nn)=ssu2(idis(nn)+1,iz(nn));
        AA3(nn)=ssu2(idis(nn),iz(nn)+1);
        AA4(nn)=ssu2(idis(nn)+1,iz(nn)+1);
      end
      sur=ps.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);

      for nn=1:nps;
        AA1(nn)=ssu3(idis(nn),iz(nn));
        AA2(nn)=ssu3(idis(nn)+1,iz(nn));
        AA3(nn)=ssu3(idis(nn),iz(nn)+1);
        AA4(nn)=ssu3(idis(nn)+1,iz(nn)+1);
      end
      sut=sh.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);

%    
%   contribution from the dip-slip comoents
      ps=pmoment(2,:)'.*co+pmoment(5,:)'.*si;
      sh=pmoment(2,:)'.*si-pmoment(5,:)'.*co;
      
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);

      for nn=1:nps;
        AA1(nn)=dsu1(idis(nn),iz(nn));
        AA2(nn)=dsu1(idis(nn)+1,iz(nn));
        AA3(nn)=dsu1(idis(nn),iz(nn)+1);
        AA4(nn)=dsu1(idis(nn)+1,iz(nn)+1);
      end
      duz=ps.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);

      for nn=1:nps;
        AA1(nn)=dsu2(idis(nn),iz(nn));
        AA2(nn)=dsu2(idis(nn)+1,iz(nn));
        AA3(nn)=dsu2(idis(nn),iz(nn)+1);
        AA4(nn)=dsu2(idis(nn)+1,iz(nn)+1);
      end
      dur=ps.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);
      for nn=1:nps;
        AA1(nn)=dsu3(idis(nn),iz(nn));
        AA2(nn)=dsu3(idis(nn)+1,iz(nn));
        AA3(nn)=dsu3(idis(nn),iz(nn)+1);
        AA4(nn)=dsu3(idis(nn)+1,iz(nn)+1);
      end
      dut=sh.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
%  contributions from the clvd component
      ps=pmoment(3,:)';
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);

      for nn=1:nps;
        AA1(nn)=clu1(idis(nn),iz(nn));
        AA2(nn)=clu1(idis(nn)+1,iz(nn));
        AA3(nn)=clu1(idis(nn),iz(nn)+1);
        AA4(nn)=clu1(idis(nn)+1,iz(nn)+1);
      end
      cluz=ps.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
      AA1=zeros(nps,1);
      AA2=zeros(nps,1);
      AA3=zeros(nps,1);
      AA4=zeros(nps,1);
      for nn=1:nps;
        AA1(nn)=clu2(idis(nn),iz(nn));
        AA2(nn)=clu2(idis(nn)+1,iz(nn));
        AA3(nn)=clu2(idis(nn),iz(nn)+1);
        AA4(nn)=clu2(idis(nn)+1,iz(nn)+1);
      end
      clur=ps.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
      
      Uz=suz+duz+cluz;
      Ur=sur+dur+clur;
      Ut=sut+dut;
    
      Ux=sum(Ur.*co-Ut.*si);
      Uy=sum(Ur.*si+Ut.*co);
      Uz=sum(Uz);
      
      Un(k)=Ux;
      Ue(k)=Uy;
      Up(k)=-Uz;
   end
  
%    G(:,ipatch+npatch)=[Ue;Un;Up];
    G(:,ipatch+npatch) = Ue.*ve + Un.*vn + Up.*vz;

end