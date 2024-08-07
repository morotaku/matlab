clear all
clc
 [z,x]=meshgrid(0:0.01:40, -0.3*10^-0:1*10^-4:0.3*10^-0); % z in cm and x in cm
 lam=1064*10^-7; % cm   wavelength
 theta=pi/180; %    radians
 n=1.5; %  refractive index
 k=2*pi/lam;
 kx=k*(n-1)*theta;
%% %%%%%%% Gaussian %%%%%%%%%%%%%%%%%%%
w0=0.1;
zR=pi.*w0.^2/lam;
G=atan(z./zR);
R=z+zR^2./z;
W=w0*sqrt(1+z.^2./zR^2);
%%
Gaus=w0^2./W.^2.*exp(-2.*(x.^2+kx.^2.*z.^2./k^2)./W.^2);
Rd=kx.*x./(1+1i.*z./zR);
 J0=besselj(0,Rd);
 J0=abs(J0);
 I0=Gaus.*J0.^2;
 I0=I0./(max(I0));
  %HGM_n=HGM./max(HGM(:));
 %mesh(z,x,I0)
 %axis off
 %I0=angle(J0);
  imshow(I0)
%  %improfile
 colormap(hot)