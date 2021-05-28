function [x_out,y_out,data_out]=multi_look(x_in,y_in,data_in,nx,ny);
% by Kang Wang in Jan. 10th, 2014
% Last Update: Kang Wang on Jan. 10th, 2014

[Ny_in,Nx_in]=size(data_in);
Nx_out = ceil(Nx_in/nx);
Ny_out = ceil(Ny_in/ny);
xout = zeros(1,Nx_out);
yout = zeros(1,Ny_out);
data_out = zeros(Ny_out,Nx_out);
for i =1:Nx_out;
   indx1= (i-1)*nx+1;
   indx2= indx1+nx;
   if (indx2>Nx_in);
     indx2=Nx_in;
   end
   x_out(i)=mean(x_in(indx1:indx2));
   for j=1:Ny_out;
      indy1=(j-1)*ny+1;
      indy2=indy1+ny;
      if (indy2>Ny_in);
        indy2=Ny_in;
      end
      y_out(j)=mean(y_in(indy1:indy2));
      data_tmp = data_in(indy1:indy2,indx1:indx2);     
      indx_good = ~isnan(data_tmp);
      data_good = data_tmp(indx_good);
      if (length(data_good)>0);
         data_out(j,i)=mean(data_good);
      else
         data_out(j,i)=NaN;
      end
   end
end

