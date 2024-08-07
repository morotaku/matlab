clear all
close all
clc
%parameters
l=3;
lam=0.63; % micrrons
k=2*pi/lam;
f=10^5; % micrrons
d=1000; %Beam waist
w0=lam*f/(pi*d) % micrrons
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
Lx=100; % micrrons
N=150;
dz=Lz/N;
M=100;
dx=Lx/M;
Z=-Lz/2:dz:Lz/2-dz;
X=-Lx/2:dx:Lx/2-dx;
[z,x]=meshgrid(Z,X);
%[z,x]=meshgrid(-3*10^-3:10^-5:3*10^-3, -1.0*10^-3:3*10^-5:1.0*10^-3); % z in mm and x in mm
%% %%%%%%%%%%%%%%%%propagation parameters %%%%%%%%%%%%%%%%%%%%
%A=-z./f;
A=1-z./f;
%B=-(z./f)*s+f+z;
B=s-(z./f)*s+z;
C=-1/f;
D=1-s/f;
A=1-z/f;
B=z;
C=-1/f;
D=1;
w=w0.*sqrt(A.^2+(B.*lam./(pi*w0^2)).^2);
Q=(C+D./q0)./(A+B./q0);  %Q=1/q
%% HG beam equation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=0;
for m=0:l
HGB=nchoosek(l,m).*(A-B./q0).^m./(A+B./q0).^(m+1).*exp(-1i.*k.*Q.*x.^2/2).*exp(-1i*k.*z).*laguerre(m,2.*x.^2./w.^2).*(-1).^m;
% HGB=nchoosek(l,m).*(A-B./q0).^m./(A+B./q0).^(m+1).*exp(-1i.*k.*Q.*x.^2/2).*exp(-1i.*k.*z).*laguerreL(m,2*x.^2./w.^2);
% H=factorial(l)/2^l.*symsum(HGB,m,0,l);
H=H+HGB;
end

HG=factorial(l)/2^l.*H;

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

function y=laguerre(n,x)
y=0;
    if n==0
        y=1;
    else
        for m=0:n
            y=y+(factorial(n)^2*(-1)^m.*x.^m)./(factorial(m)*factorial(m)*factorial(n-m));
        end
    end
end