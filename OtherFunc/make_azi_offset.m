%doing the cross-correlation of two amplitudes images
%by Kang Wang in May 2016
clc
clear
format long

PRMFILE='xcorr_matlab.PRM';
% 
grd_file1=load_PRM(PRMFILE,'GRD_FILE1');
grd_file2=load_PRM(PRMFILE,'GRD_FILE2');

disp(['Reading data start ....']);
[x1,y1,z1]=grdread2(grd_file1);
[x2,y2,z2]=grdread2(grd_file2);

xxmin=min(x1)-1;
yymin=min(y1)-1;


z1=z1*1.0e5;
z2=z2*1.0e5;

disp(['Reading data finish ....']);
[ny,nx]=size(z1);

xmin=1;
xmax=nx;

ymin=1;
ymax=ny;

Nx_output=load_PRM(PRMFILE,'NX_OUTPUT');
Ny_output=load_PRM(PRMFILE,'NY_OUTPUT');
nx_search=load_PRM(PRMFILE,'NX_SEARCH');
ny_search=load_PRM(PRMFILE,'NY_SEARCH');

xout_min=xmin+floor(nx_search/2);
xout_max=xmax-floor(ny_search/2);
xout=floor(linspace(xout_min,xout_max,Nx_output));

yout_min=ymin+floor(ny_search/2);
yout_max=ymax-floor(ny_search/2);
yout=floor(linspace(yout_min,yout_max,Ny_output));

f_xcorr=fopen('xcorr.dat','w');
for i=1:Ny_output;
    for j=1:Nx_output;
        xout_this_box=xout(j);
        yout_this_box=yout(i);
        indx_y_min=yout_this_box-ny_search/2;
        indx_y_max=yout_this_box+ny_search/2-1;
        
        indx_x_min=xout_this_box-nx_search/2;
        indx_x_max=xout_this_box+nx_search/2-1;
        
        img1=z1(indx_y_min:indx_y_max,indx_x_min:indx_x_max);
        img2=z2(indx_y_min:indx_y_max,indx_x_min:indx_x_max);
        m1=mean(mean(img1));
        m2=mean(mean(img2));
        img1_new = img1-m1;
        img2_new = img2-m2;
%        [output Greg] = dftregistration(fft2(img2_new),fft2(img1_new),1000);
        [output Greg] = xcorr_insar(fft2(img2_new),fft2(img1_new),100);
       
        dx=output(4);
        dy=output(3);
        snr=output(1);
%        fprintf(f_xcorr,'%10d %8.3f %10d %8.3f %8.2f\n',[xout_this_box,dx,yout_this_box,dy,snr]);
        fprintf(f_xcorr,'%12.2f %8.3f %12.2f %8.3f %8.2f\n',[xout_this_box+xxmin,dx,yout_this_box+yymin,dy,snr]);
    end
       disp(['Done Percentage(%) : ',num2str(100*i/Ny_output)]);
end
fclose(f_xcorr);
quit;
