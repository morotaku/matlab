close all
clear all
clc
l=0;
w0=0.01; % microns
z=0.1;
lam=1064*10^-7; % cm   wavelength
theta=1*pi/180; %    radians
n=1.5; %  refractive index
k=2*pi/lam;
kr=k*(n-1)*theta;
%% %%%%%%%%%%%%%$$$$$$$$$$$$$$%%%%%%%%%%%%
[x, y]=meshgrid(-3*10^-1:10^-3:3*10^-1, -3*10^-1:10^-3:3*10^-1);
r=x.^2+y.^2;
theta=atan(y./x);
R=kr.*r;
A=sqrt(k.^2.*r.^2./z.^2);
B=k.*w0.^2./(k.*w0.^2+2.*1i.*z);
C=1i.*kr.*A.*w0.^2.*z./(k.*w0.^2+2.*1i.*z);
D=(kr.^2+A.^2).*w0.^2./(2.*k.^2.*w0.^2+4.*1i.*k.*z);
%% 
U1=exp(1i.*k.*(z-kr.^2./(2.*k.^2))+1i.*l.*pi).*besselj(l,R);
U2=B.*besseli(l,C).*exp(1i.*l.*pi/2).*exp(1i.*k.*(1+r.^2./(2.*z.^2)-D));
I=abs(U1-U2).^2;
 I0=I/max(max(I));
 imshow(I0)