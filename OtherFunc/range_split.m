clc
clear

prmfile='range_split.PRM';
GMTSAR_PRM=load_PRM(prmfile,'GMTSAR_PRM')


slc_file=load_gmtsar_PRM(GMTSAR_PRM,'SLC_file');
nx=load_gmtsar_PRM(GMTSAR_PRM,'num_rng_bins');
ny=load_gmtsar_PRM(GMTSAR_PRM,'nrows');

ndata_line = 2*nx;

filt_low=ones(nx,1);
filt_high=ones(nx,1);

filt_low(1:floor(nx/2))=0;
filt_high(floor(nx/2):nx)=0;

slc_h=[slc_file,'H'];
slc_l=[slc_file,'L'];

fc=fopen(slc_file,'r');
fh=fopen(slc_h,'wb');
fl=fopen(slc_l,'wb');

slc_row_h=zeros(ndata_line,1);
slc_row_l=zeros(ndata_line,1);

for i=1:ny;
   row_readin=fread(fc,[ndata_line,1],'int16');
   slc_real=row_readin(1:2:ndata_line);
   slc_imag=row_readin(2:2:ndata_line);
   slc_row=complex(slc_real,slc_imag);
  
   slc_spec=fft(slc_row);
   slc_h=ifft(slc_spec.*filt_high);
   slc_l=ifft(slc_spec.*filt_low);

   rh=real(slc_h);
   ih=imag(slc_h);

   rl=real(slc_l);
   il=imag(slc_l);
   
   slc_row_h(1:2:ndata_line)=rh;
   slc_row_h(2:2:ndata_line)=ih;
  
   slc_row_l(1:2:ndata_line)=rl;
   slc_row_l(2:2:ndata_line)=il;
  
   fwrite(fh,slc_row_h,'int16');
   fwrite(fl,slc_row_l,'int16');
   if (mod(i,500)==0);
     disp(['Percentage Done (%): ',num2str(100*i/ny)]);
   end
end

fclose(fc);
fclose(fh);
fclose(fl);

quit;
