function G=calc_green_lines_edcmp_optimize2(slip_model,data_insar,finp)
% Calculate the Green's function for InSAR observations using EDGRN
% by Kang Wang on 08/21/2019
% Last Updated by Kang Wang on 02/20/2020
% slip_patch is used to scale the line sources moment

format long

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
%superposition of all discrete point sources
Ue_ss=zeros(nobs,1);
Ue_ds=zeros(nobs,1);

Un_ss=zeros(nobs,1);
Un_ds=zeros(nobs,1);

Up_ss=zeros(nobs,1);
Up_ds=zeros(nobs,1);

ZSEPS=1.0d-2;

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
% clu3=zeros(nr+2,nz+2);
clu1(1:nr,1:nz)=reshape(cl.uz,nr,nz);
clu2(1:nr,1:nz)=reshape(cl.ur,nr,nz);
% clu3(1:nr,1:nz)=reshape(cl.ut,nr,nz);


dr=(r2-r1)/(nr-1);
dz=(z2-z1)/(nz-1);

if ((zrec0>=z1) && (zrec0<=z2) && (r1==0))
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

    for ipatch = 1:npatch
       if mod(ipatch, 400) == 0
           disp(['... Green-GPS patch: ',' [',num2str(ipatch),'/',num2str(npatch),']']);
       end

       xs_patch=slip_model(ipatch,2);
       ys_patch=slip_model(ipatch,1);
       zs_patch=-slip_model(ipatch,3);  % reverse the z-axis
       lp_patch=slip_model(ipatch,4);
       wp_patch=slip_model(ipatch,5);
       strike_patch=slip_model(ipatch,6);
       dip_patch=slip_model(ipatch,7);
       slip_patch=slip_model(ipatch,8);  % scaled by the moment

       [pxs_ss,pys_ss,pzs_ss,pmoment_ss]=edcdisc(xs_patch,ys_patch,zs_patch,slip_patch,lp_patch,...
       wp_patch,strike_patch,dip_patch,0,nz,z1,z2,dr,dz);

       [~, ~, ~, pmoment_ds]=edcdisc(xs_patch,ys_patch,zs_patch,slip_patch,lp_patch,...
       wp_patch,strike_patch,dip_patch,90,nz,z1,z2,dr,dz);   

       nps=length(pxs_ss);   %number of point sources
       AA1=zeros(nps,1);
       AA2=zeros(nps,1);
       AA3=zeros(nps,1);
       AA4=zeros(nps,1);

       for k=1:nobs
          dis=sqrt((xrec(k)-pxs_ss).^2+(yrec(k)-pys_ss).^2);
          if (dis<0)
              azi=0;
          else
              azi=atan2(yrec(k)-pys_ss,xrec(k)-pxs_ss);  
          end

          iz=floor((pzs_ss-z1)/dz)+1;
          dzs=(pzs_ss-(z1+dz*(iz-1)))/dz;

          if (dis<r1)
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

    %  contributions from the strike-slip components
          ps_ss=pmoment_ss(1,:)'.*si2+pmoment_ss(4,:)'.*co2;
          sh_ss=pmoment_ss(1,:)'.*co2-pmoment_ss(4,:)'.*si2;
          ps_ds=pmoment_ds(1,:)'.*si2+pmoment_ds(4,:)'.*co2;
          sh_ds=pmoment_ds(1,:)'.*co2-pmoment_ds(4,:)'.*si2;

          for nn=1:nps
            AA1(nn)=ssu1(idis(nn),iz(nn));
            AA2(nn)=ssu1(idis(nn)+1,iz(nn));
            AA3(nn)=ssu1(idis(nn),iz(nn)+1);
            AA4(nn)=ssu1(idis(nn)+1,iz(nn)+1);
          end
          suz_ss=ps_ss.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
          suz_ds=ps_ds.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);

          for nn=1:nps
            AA1(nn)=ssu2(idis(nn),iz(nn));
            AA2(nn)=ssu2(idis(nn)+1,iz(nn));
            AA3(nn)=ssu2(idis(nn),iz(nn)+1);
            AA4(nn)=ssu2(idis(nn)+1,iz(nn)+1);
          end     
          sur_ss=ps_ss.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
          sur_ds=ps_ds.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);

          for nn=1:nps
            AA1(nn)=ssu3(idis(nn),iz(nn));
            AA2(nn)=ssu3(idis(nn)+1,iz(nn));
            AA3(nn)=ssu3(idis(nn),iz(nn)+1);
            AA4(nn)=ssu3(idis(nn)+1,iz(nn)+1);
          end    
          sut_ss=sh_ss.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
          sut_ds=sh_ds.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);

    %   contribution from the dip-slip comoents
          ps_ss=pmoment_ss(2,:)'.*co+pmoment_ss(5,:)'.*si;
          sh_ss=pmoment_ss(2,:)'.*si-pmoment_ss(5,:)'.*co;
          ps_ds=pmoment_ds(2,:)'.*co+pmoment_ds(5,:)'.*si;
          sh_ds=pmoment_ds(2,:)'.*si-pmoment_ds(5,:)'.*co;

          for nn=1:nps
            AA1(nn)=dsu1(idis(nn),iz(nn));
            AA2(nn)=dsu1(idis(nn)+1,iz(nn));
            AA3(nn)=dsu1(idis(nn),iz(nn)+1);
            AA4(nn)=dsu1(idis(nn)+1,iz(nn)+1);
          end
          duz_ss=ps_ss.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
          duz_ds=ps_ds.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);

          for nn=1:nps
            AA1(nn)=dsu2(idis(nn),iz(nn));
            AA2(nn)=dsu2(idis(nn)+1,iz(nn));
            AA3(nn)=dsu2(idis(nn),iz(nn)+1);
            AA4(nn)=dsu2(idis(nn)+1,iz(nn)+1);
          end
          dur_ss=ps_ss.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
          dur_ds=ps_ds.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);

          for nn=1:nps
            AA1(nn)=dsu3(idis(nn),iz(nn));
            AA2(nn)=dsu3(idis(nn)+1,iz(nn));
            AA3(nn)=dsu3(idis(nn),iz(nn)+1);
            AA4(nn)=dsu3(idis(nn)+1,iz(nn)+1);
          end
          dut_ss=sh_ss.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
          dut_ds=sh_ds.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);

    %  contributions from the clvd component
          ps_ss=pmoment_ss(3,:)';
          ps_ds=pmoment_ds(3,:)';

          for nn=1:nps
            AA1(nn)=clu1(idis(nn),iz(nn));
            AA2(nn)=clu1(idis(nn)+1,iz(nn));
            AA3(nn)=clu1(idis(nn),iz(nn)+1);
            AA4(nn)=clu1(idis(nn)+1,iz(nn)+1);
          end
          cluz_ss=ps_ss.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
          cluz_ds=ps_ds.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);

          for nn=1:nps
            AA1(nn)=clu2(idis(nn),iz(nn));
            AA2(nn)=clu2(idis(nn)+1,iz(nn));
            AA3(nn)=clu2(idis(nn),iz(nn)+1);
            AA4(nn)=clu2(idis(nn)+1,iz(nn)+1);
          end
          clur_ss=ps_ss.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);
          clur_ds=ps_ds.*(w00.*AA1+w10.*AA2+w01.*AA3+w11.*AA4);

          Un_ss(k) = sum((sur_ss+dur_ss+clur_ss).*co - (sut_ss+dut_ss).*si);
          Ue_ss(k) = sum((sur_ss+dur_ss+clur_ss).*si + (sut_ss+dut_ss).*co);
          Up_ss(k) = -sum(suz_ss+duz_ss+cluz_ss);

          Un_ds(k) = sum((sur_ds+dur_ds+clur_ds).*co - (sut_ds+dut_ds).*si);
          Ue_ds(k) = sum((sur_ds+dur_ds+clur_ds).*si + (sut_ds+dut_ds).*co);
          Up_ds(k) = -sum(suz_ds+duz_ds+cluz_ds);
       end

       G(:,ipatch) = Ue_ss.*ve + Un_ss.*vn + Up_ss.*vz;
       G(:,ipatch+npatch) = Ue_ds.*ve + Un_ds.*vn + Up_ds.*vz;
    end
  
end