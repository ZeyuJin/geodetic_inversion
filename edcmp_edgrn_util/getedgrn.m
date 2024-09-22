function getedgrn(finp);
% by Kang Wang on 02/20/20;
% finp is the EDGRN input file;

fid = fopen(finp,'r');
igood=0;
iline=0;
while ~feof(fid)
A = fgetl(fid);
 if (~ismember('#',A))
  iline=iline+1;
  if (iline==5);
    B=split(A,"'");
    dir_grn=B{2};
    fss=B{4};
    fds=B{6};
    fcl=B{8};
   break;
  end
 end
end
fclose(fid);

fss=[dir_grn,fss];
fds=[dir_grn,fds];
fcl=[dir_grn,fcl];

fid=fopen(fss,'r');
ncomment=0;
while ~feof(fid)
  A=fgetl(fid);
  ncomment=ncomment+1;
  if (~ismember('#',A));
   A=str2num(A);
   nr=A(1);
   r1=A(2);
   r2=A(3);
   nz=A(4);
   z1=A(5);
   z2=A(6);
   zrec0=A(7);
   lambda=A(8);
   mu=A(9);
      break;
  end
end
fclose(fid);

[uz,ur,ut,ezz,err,ett,ezr,ert,etz,duz_dr]=textread(fss,'%f %f %f %f %f %f %f %f %f %f\n ','headerlines',ncomment);
%save('edgrn_ss.mat','nr','nz','r1','r2','z1','z2','zrec0','lambda','mu',...
%    'uz','ur','ut','ezz','err','ett','ezr','ert','etz','duz_dr');

save(['edgrn_ss_',finp,'.mat'],'nr','nz','r1','r2','z1','z2','zrec0','lambda','mu',...
    'uz','ur','ut','ezz','err','ett','ezr','ert','etz','duz_dr');

fid=fopen(fds,'r');
ncomment=0;
while ~feof(fid);
  A=fgetl(fid);
  ncomment=ncomment+1;
  if (~ismember('#',A));
   A=str2num(A);
   nr=A(1);
   r1=A(2);
   r2=A(3);
   nz=A(4);
   z1=A(5);
   z2=A(6);
   zrec0=A(7);
   lambda=A(8);
   mu=A(9);
      break;
  end
end
fclose(fid);

[uz,ur,ut,ezz,err,ett,ezr,ert,etz,duz_dr]=textread(fds,'%f %f %f %f %f %f %f %f %f %f\n ','headerlines',ncomment);
%save('edgrn_ds.mat','nr','nz','r1','r2','z1','z2','zrec0','lambda','mu',...
%    'uz','ur','ut','ezz','err','ett','ezr','ert','etz','duz_dr');


save(['edgrn_ds_',finp,'.mat'],'nr','nz','r1','r2','z1','z2','zrec0','lambda','mu',...
    'uz','ur','ut','ezz','err','ett','ezr','ert','etz','duz_dr');


fid=fopen(fcl,'r');
ncomment=0;
while ~feof(fid);
  A=fgetl(fid);
  ncomment=ncomment+1;
  if (~ismember('#',A));
   A=str2num(A);
   nr=A(1);
   r1=A(2);
   r2=A(3);
   nz=A(4);
   z1=A(5);
   z2=A(6);
   zrec0=A(7);
   lambda=A(8);
   mu=A(9);
      break;
  end
end
fclose(fid);

[uz,ur,ut,ezz,err,ett,ezr,ert,etz,duz_dr]=textread(fcl,'%f %f %f %f %f %f %f %f %f %f\n ','headerlines',ncomment);
%save('edgrn_cl.mat','nr','nz','r1','r2','z1','z2','zrec0','lambda','mu',...
%    'uz','ur','ut','ezz','err','ett','ezr','ert','etz','duz_dr');

save(['edgrn_cl_',finp,'.mat'],'nr','nz','r1','r2','z1','z2','zrec0','lambda','mu',...
    'uz','ur','ut','ezz','err','ett','ezr','ert','etz','duz_dr');
