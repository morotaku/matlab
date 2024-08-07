clear all
close all
clc
format shortEng

%%% x-y　coordinate
N=200;
Lx=300; %Display lange
Lz=15000; 
X=linspace(-Lx,Lx,N);
Z=linspace(0,15000,N);
[z,x]=meshgrid(Z,X);

%%% Topological charge
l=3;
%%% Radial index
n=2;

%%% Parameters(µm)
d=1000; %Incident beam size
lam=0.64; %Wave length
k=2*pi/lam; %Wave number
f=10^5; %Focus length
w0=1000;%Beam waist
w=lam*f/(pi*d)
zr=w^2*k/2 %Rayleigh length
s=0; %Distance from the input plane to lens

G0=1;
i=1j;

%ABCD matrix
A=-z./f;
B=-(z./f).*s+f+z;
C=-1/f;
D=1-s/f;
R=z+zr^2./z; %Beam curvature
q0=1j*k*w0^2/2;
Q=(A+B./q0)./(C+D./q0);
W=w0.*sqrt(A.^2+(B.*lam./(pi*w0^2)).^2);
%%%HGB
C1=1j*k*G0*factorial(n)./(2*w0^(2*n).*B).*((1./w0.^2)+(1j*k.*A./(2.*B))).^(-1-n).*exp(-1j*k.*z);
C2=exp(-1j*k.*(D.*x.^2)./(2.*B));
C3=exp(-(k.*x./(2.*B)).^2./((1./w0.^2)+(1j*k.*A./(2.*B)))).*laguerre(n,(k.*x./(2.*B)).^2./((1./w0.^2)+(1j*k.*A./(2.*B))));
hgb=C1.*C2.*C3;
HGB=abs(hgb.^2);
HGB_n=HGB./max(HGB(:));

imagesc(Z,X,HGB_n,[0,1])
colormap("hot")
%xlim([0 4000]);
%shading interp; lighting phong; view(2); axis equal; axis tight;
%axis off;

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