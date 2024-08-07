clear all
close all
clc
format shortEng

%%% x-y　coordinate
N=401;
L=4000; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);
Z=linspace(-6000,2000,N);

%%% Topological charge
l=0;
%%% Radial index
n=1;


%%% Parameters(µm)
w=1; %Beam waist
lam=1.06; %Wave length
k=2*pi/lam; %Wave number
zr=w^2*k/2 %Rayleigh length
f=1000; %Focus length
s=0; %Distance from the input plane to lens

G0=1;
i=1j;

for j=0:400
    z=20*j-6000; %Beam position
    %ABCD matrix
    A=1-z/f;
    B=s+z-z*s/f;
    %A=-z/f;
    %B=-(z/f)*s+f+z;
    C=-1/f;
    D=1-s/f; 
    R=z+zr^2/z; %Beam curvature
    q0=1j*k*w^2/2;
    Q=(A+B/q0)/(C+D/q0);
    W=w*sqrt(A^2+(B*lam/(pi*w^2))^2);
    %%%HGB
    C1=G0*factorial(n)/2^n;
    C2=zeros(N,N);
    for m=0:n
        C2=C2+(-1)^m*nchoosek(n,m)*((A-B/q0)^m/(A+B/q0)^(m+1))*exp(-1j*k*z).*laguerre(m,2/w^2.*r.^2).*exp(-1j*k/(2*Q).*r.^2);
    end
    lag=C1.*C2;
    HGM=abs(lag.^2);
    HGM_n=HGM./max(HGM(:));
    %HGM_n(:,500)
    HGM_xz(:,j+1)=HGM_n(:,200);
end
    imagesc(Z,X,HGM_xz)
    colormap("hot")
    %xlim([0 4000]);
    %shading interp; lighting phong; view(2); axis equal; axis tight; 
    %axis off;
