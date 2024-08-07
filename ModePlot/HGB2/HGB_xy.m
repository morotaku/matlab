clear all
close all
clc
format shortEng
z=input('z=')
%%% x-y　coordinate
N=200;
L=20; 
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);
%%% Topological charge
l=3;
%%% Radial index
n=3;


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
A=-z/f;
B=-(z/f)*s+f+z;
C=-1/f;
D=1-s/f;
R=z+zr^2./z; %Beam curvature
q0=1j*k*w0^2/2;
Q=(A+B./q0)./(C+D./q0);
W=w0.*sqrt(A.^2+(B.*lam./(pi*w0^2)).^2)
%%%HGB
C1=1j*k*G0*factorial(n)./(2*w0^(2*n).*B).*((1./w0.^2)+(1j*k.*A./(2.*B))).^(-1-n).*exp(-1j*k.*z);
C2=exp(-1j*k.*(D.*r.^2)./(2.*B));
C3=exp(-(k.*r./(2.*B)).^2./((1./w0.^2)+(1j*k.*A./(2.*B)))).*laguerre(n,(k.*r./(2.*B)).^2./((1./w0.^2)+(1j*k.*A./(2.*B))));
hgb=C1.*C2.*C3;
hgb_n=abs(hgb)./max(max(abs(hgb)));
HGB=hgb_n.^2;
HGB_n=HGB./max(HGB(:));

figure(1)
imagesc(X,Y,HGB_n,[0,1])
colormap("hot")
axis equal; axis tight;
figure(2)
plot(X,HGB_n(:,N/2))
%xlim([0 4000]);
%shading interp; lighting phong; view(2); 
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