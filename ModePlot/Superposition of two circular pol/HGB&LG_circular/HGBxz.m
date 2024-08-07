clear all
close all
clc
format shortEng

%%% x-y　coordinate
N=1001;
L=1*10^2; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);
Z=linspace(-1*10^3,6*10^3,N);

%%% Topological charge
l=0;
%%% Radial index
n=3;


%%% Parameters(µm)
w=1000; %Beam waist
lam=0.63; %Wave length
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
f=10^5; %Focus length
%s=10^6; %Distance from the input plane to lens
w2=lam*f/(pi*w)
zr2=pi*w2^2/lam
G0=1;
i=1j;

for j=0:300
    z=7*(j/3)*10^1-10^3; %Beam position
    %ABCD matrix
    A=1-z/f;
    B=z;
    C=-1/f;
    D=1;
    %A=1;
    %B=z;
    %C=0;
    %D=1;
    R=z+zr^2/z; %Beam curvature
    
    q0=1j*k*w2^2/2;
    Q=(A+B/q0)/(C+D/q0);
    zr2=pi*w2^2/lam;
    W=w2*sqrt(A^2+(B*lam/(pi*w2^2))^2);
    %%%Mode plot
    C1=G0*factorial(n)/2^n;
    C2=0;
 
    for m=0:n
        C2=C2+(-1)^m*nchoosek(n,m)*((A-B/q0)^m/(A+B/q0)^(m+1))*exp(-1j*k*z).*laguerre(m,(2.*r.^2)./W^2).*exp(-1j*k/(2*Q).*r.^2);
        %C2=C2+(-1)^m*nchoosek(n,m)*((A-B/q0)^m/(A+B/q0)^(m+1)).*laguerre(m,2/w^2.*r.^2).*exp(-1j*k/(2*Q).*r.^2);
    end
    hgb=C1.*C2;
    HGM=abs(hgb.^2);
    HGM_n=HGM./max(HGM(:));
    %HGM_n(:,500)
    HGM_xz(:,j+1)=HGM_n(:,500);
end
    imagesc(Z,X,HGM_xz)
    colormap("hot")
    %xlim([0 4000]);
    %shading interp; lighting phong; view(2); axis equal; axis tight; 
    %axis off;
