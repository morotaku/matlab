clear all
close all
clc
format shortEng

%%% x-y　coordinate
N=401;
L=20; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);
Z=linspace(0,40,N);

%%% Topological charge
l=0;
%%% Radial index
p=2;


%%% Parameters(µm)
w=10000; %Beam waist
lam=1.06; %Wave length
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
f=20000; %Focus length
s=200000; %Distance from the input plane to lens
G0=1;
i=1j;

for j=0:400
    z=j/10; %Beam position
    R=z+zr^2/z; %Beam curvature
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    
    %ABCD matrix
    A=-z/f;
    B=-(z/f)*s+f+z;
    C=-1/f;
    D=1-s/f; 
    
    %%%HGB
    C1=(i*k*G0*factorial(p)/(2*B*w^(2*p)))*((1/(w^2)+(i*k*A/(2*B)))^(-1-p))*exp(-i*k*z);
    C2=exp(((-i*k*D.*(r.^2))./(2*B)));
    C3=exp(-((k.*r)./(2*B)).^2./(1/w^2+(i*k*A/(2*B))));
    lag=laguerre(p,(k/(2*B))^2/((1/w^2)+(i*k*A/(2*B))).*(r.^2));
    hgb=C1.*C2.*C3.*lag;
    HGM=abs(hgb.^2);
    HGM_n=HGM./max(HGM(:));
    %HGM_n(:,500)
    HGM_xz(:,j+1)=HGM_n(:,200);
end
    imagesc(Z,X,HGM_xz)
    colormap("hot")
    shading interp; lighting phong; view(2); axis equal; axis tight; 
    %axis off;
