%Last Updated by Kang Wang on 12/21/2016
clc
clear

SOL=299792458; %speed of light

[prm_list]=textread('SLC.list','%s\n');
Nprm=length(prm_list);
for nn=1:Nprm;
    this_prm=prm_list{nn};
    disp(['working on ',this_prm]);
    GMTSAR_PRM=[this_prm,'.PRM'];
    
    slc_file=load_gmtsar_PRM(GMTSAR_PRM,'SLC_file');
    nx=load_gmtsar_PRM(GMTSAR_PRM,'num_rng_bins');
    ny=load_gmtsar_PRM(GMTSAR_PRM,'nrows');
    nr=load_gmtsar_PRM(GMTSAR_PRM,'near_range');
    rng_samp_rate=load_gmtsar_PRM(GMTSAR_PRM,'rng_samp_rate');
    wavelength=load_gmtsar_PRM(GMTSAR_PRM,'radar_wavelength');
    f0=SOL/wavelength;  %center frequency


    ndata_line = 2*nx;
    filter_H=zeros(nx,1);
    filter_L=zeros(nx,1);
    filter_H(1:floor(nx/2))=1;
    filter_L(floor(nx/2):nx)=1;


    slch_file=[slc_file,'H'];
    slcl_file=[slc_file,'L'];

    f_slc=fopen(slc_file,'r');
    f_slc_h=fopen(slch_file,'wb');
    f_slc_l=fopen(slcl_file,'wb');

    slc_row_h=zeros(ndata_line,1);
    slc_row_l=zeros(ndata_line,1);

    disp(['Input SLC is: ', slc_file]);
    disp(['Output SLCH is: ', slch_file]);
    disp(['Output SLCL is: ', slcl_file]);

    Nmax=ny;
   for k=1:Nmax;
%   disp(['working on line: ', num2str(k),'/',num2str(ny)]);   
     row_readin=fread(f_slc,[ndata_line,1],'int16');
     slc_real=row_readin(1:2:ndata_line);
     slc_imag=row_readin(2:2:ndata_line);
     slc_row=complex(slc_real,slc_imag);

     spec=fft(slc_row);
     spec_L=spec.*filter_L;
     spec_H=spec.*filter_H;
  

     s1=fftshift(spec);
     s2=fftshift(spec_L);
     s3=fftshift(spec_H);
   
     A1=abs(s1).^2;
     A2=abs(s2).^2;
     A3=abs(s3).^2;
   
     indx_good1=find(A1>0.03*max(A1));
     indx_good2=find(A2>0.03*max(A2));
     indx_good3=find(A3>0.03*max(A3));
     x1_center=min(indx_good1);
     x2_center=max(indx_good1);
     xc=mean([x1_center,x2_center]);
   
     x1_high=min(indx_good3);
     x2_high=max(indx_good3);
     xc_high=mean([x1_high,x2_high]);
   
     x1_low=min(indx_good2);
     x2_low=max(indx_good2);
     xc_low=mean([x1_low,x2_low]);
   
     Nh_shift=fix(xc_high-xc);
     Nl_shift=fix(xc_low-xc);
   
     spec_L_new=circshift(spec_L,Nh_shift);
     spec_H_new=circshift(spec_H,Nl_shift); 
     slc_H=ifft(spec_H_new);
     slc_L=ifft(spec_L_new);  

     rh=real(slc_H); %real part of High frequency
     ih=imag(slc_H); %imaginary part of High frequency
    
     rl=real(slc_L);
     il=imag(slc_L);

     slc_row_h(1:2:ndata_line)=rh;
     slc_row_h(2:2:ndata_line)=ih;

     slc_row_l(1:2:ndata_line)=rl;
     slc_row_l(2:2:ndata_line)=il;

     fwrite(f_slc_h,slc_row_h,'int16');
     fwrite(f_slc_l,slc_row_l,'int16');

    if (mod(k,500)==0);
       disp(['Percentage Done (%): ',num2str(100*k/ny)])
    end

    if (k==Nmax);
       disp(['Percentage Done (%): ',num2str(100*k/Nmax)])
    end

   end

 fclose(f_slc);
 fclose(f_slc_l);
 fclose(f_slc_h);
end
quit;
