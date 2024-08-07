clear all
close all
clc
%syms m
%parameters
l=1;
lam=1065*10^-3; % micrrons
k=2*pi/lam;
f=5*10^3; % micrrons
w0=572*10^-3; % micrrons
q0=1i*k*w0^2/2;
s=100;
% lam=1065*10^-6; % mm
%  k=2*pi/lam;
%  f=5; % mm
%  w0=450*10^-6; % mm
% q0=1i*k*w0^2/2;
%% %%%%%%%% coordinates %%%%%%%%%%%%
Lz=5000; % micrrons
Lx=5000; % micrrons
N=200;
dz=Lz/N;
M=100;
dx=Lx/M;
Z=-Lz/2:dz:Lz/2-dz;
X=-Lx/2:dx:Lx/2-dx;
[z,x]=meshgrid(Z,X);
%[z,x]=meshgrid(-3*10^-3:10^-5:3*10^-3, -1.0*10^-3:3*10^-5:1.0*10^-3); % z in mm and x in mm
%% %%%%%%%%%%%%%%%%propagation parameters %%%%%%%%%%%%%%%%%%%%
A=1-z./f;
B=-(z./f)*s+f+z;
C=-1/f;
D=1-s/f;
w=w0.*sqrt(A.^2+(B.*lam./(pi*w0^2)).^2);
Q=(C+D./q0)./(A+B./q0);  %Q=1/q
%% HG beam equation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HGB=0;
for m=0:l
HGB=HGB+nchoosek(l,m).*(A-B./q0).^m./(A+B./q0).^(m+1).*exp(-1i.*k.*Q.*x.^2/2).*exp(-1i.*k.*z).*laguerre(2,2.*x.^2./w.^2).*(-1).^m;
% HGB=nchoosek(l,m).*(A-B./q0).^m./(A+B./q0).^(m+1).*exp(-1i.*k.*Q.*x.^2/2).*exp(-1i.*k.*z).*laguerreL(m,2*x.^2./w.^2);
% H=factorial(l)/2^l.*symsum(HGB,m,0,l);
%H=factorial(l)/2^l.*HGB;
end
H=factorial(l)/2^l;
HG=HGB.*H;
I=abs(HG).^2;
I=I/max(max(I));
%I=I.^2;
%A=angle(HG);
 imagesc(Z,X,I)
 colormap hot
 xlabel('z (\mum)'); ylabel('x (\mum)');
 %improfile
% %imagesc(A)
% axis equal
% axis off