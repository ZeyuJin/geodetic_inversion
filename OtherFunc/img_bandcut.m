function img_out=img_bandcut(img_in,zmin,zmax);
%by Kang Wang on 02/02/2014
%
% zmin : minimum of energy window of interest
% zmax : maximum of energy window of interest

%%%%% FFT the original image
z_fft=fftshift(fft2(img_in));

% Get the energy of the spectrum
rspe=abs(z_fft);

flag1=rspe>zmin;
flag2=rspe<zmax;

% Determine the window in the freqency domain
F_window = flag1.*flag2;

% Set the energy outside of the window to be zero
% Only consider the amplitude
my_filter=complex(F_window,ones(size(F_window)));
spe_out = z_fft.*my_filter;
img_out=real(ifft2(ifftshift(spe_out)));



