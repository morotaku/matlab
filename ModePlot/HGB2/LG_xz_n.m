clear all
close all
clc
%parameters
p=0;
l=0;
lam=0.63; % micrrons
k=2*pi/lam;
f=10^5; % micrrons
d=1000; %Beam waist
w0=20.05 % micrrons
zr=w0^2*k/2
q0=1i*k*w0^2/2;
s=0;
% lam=1065*10^-6; % mm
%  k=2*pi/lam;
%  f=5; % mm
%  w0=450*10^-6; % mm
% q0=1i*k*w0^2/2;
%% %%%%%%%% coordinates %%%%%%%%%%%%
Lz=10000; % micrrons
Lx=600; % micrrons
N=150;
dz=Lz/N;
M=100;
dx=Lx/M;
Z=0:dz:Lz;
X=-Lx/2:dx:Lx/2-dx;
[z,x]=meshgrid(Z,X);
%[z,x]=meshgrid(-3*10^-3:10^-5:3*10^-3, -1.0*10^-3:3*10^-5:1.0*10^-3); % z in mm and x in mm
%% %%%%%%%%%%%%%%%%propagation parameters %%%%%%%%%%%%%%%%%%%%
w=w0.*(1+(z./zr).^2).^0.5; %Beam size
%% LG beam equation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
%Gouy=-1j.*((-(l*phi)-((2*p+abs(l)+1)*atan(z/zr)))+k*(r.^2/(2*R))+k*z);
lg=C./w.*exp((-1).*x.^2./w.^2).*(2^0.5.*x./w).^abs(l).*laguerre(p,l,2.*x.^2./w.^2);%.*exp(Gouy);


LG=abs(lg.^2);
I=LG/max(max(LG));
%I=I.^2;
%A=angle(HG);
 imagesc(Z,X,I)
 colormap hot
 xlabel('z (\mum)'); ylabel('x (\mum)');
 %improfile
% %imagesc(A)
% axis equal



function y=laguerre(p,l,x)
y=zeros(p+1,1);
if p==0
    y=1;
else
for m=0:p
    y(p+1-m)=((-1).^m.*(factorial(p+l)))./(factorial(p-m).*factorial(l+m).*factorial(m));
end
end
y=polyval(y,x);
end