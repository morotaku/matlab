clear all
close all
clc
format shortEng

%%% x-y　coordinate
N=401;
L=300; %Display lange
Lz=10000;
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);
Z=linspace(0,Lz,N);

%%% Topological charge
l=0;
%%% Radial index
n=10;


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

for j=0:400
    z=j/400*Lz; %Beam position
    R=z+zr^2/z; %Beam curvature
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    
    %ABCD matrix
    A=-z./f;
    B=-(z./f).*s+f+z;
    C=-1/f;
    D=1-s/f; 
    
    %%%HGB
    C1=1j*k*G0*factorial(n)./(2*w0^(2*n).*B).*((1./w0.^2)+(1j*k.*A./(2.*B))).^(-1-n).*exp(-1j*k.*z);
    C2=exp(-1j*k.*(D.*r.^2)./(2.*B));
    C3=exp(-(k.*r./(2.*B)).^2./((1./w0.^2)+(1j*k.*A./(2.*B)))).*laguerre(n,(k.*r./(2.*B)).^2./((1./w0.^2)+(1j*k.*A./(2.*B))));
    hgb=C1.*C2.*C3;
    HGM=abs(hgb.^2);
    HGM_n=HGM./max(HGM(:));
    %HGM_n(:,500)
    HGM_xz(:,j+1)=HGM_n(:,200);
end
    imagesc(Z,X,HGM_xz)
    colormap("hot")
    shading interp; lighting phong; view(2); axis tight;  
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