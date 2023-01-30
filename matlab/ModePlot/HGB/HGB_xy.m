clear all
close all
clc
format shortEng
%%% Topological charge
l=0;
%%% Radial index
p=1;

%%% Parameters(µm)
w=10000; %Beam waist
lam=1.06; %Wave length
z=input('position: '); %Beam position
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size

f=20000; %Focus length
s=200000; %Distance from the input plane to lens
G0=1;
i=1j;



%ABCD matrix
A=-z/f;
B=-(z/f)*s+f+z;
C=-1/f;
D=1-s/f;

%%% x-y　coordinate
N=1001;
L=20; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);



%%%HGB
C1=(i*k*G0*factorial(p)/(2*B*w^(2*p)))*((1/(w^2)+(i*k*A/(2*B)))^(-1-p))*exp(-i*k*z);
C2=exp(((-i*k*D.*(r.^2))./(2*B)));
C3=exp(-((k.*r)./(2*B)).^2./(1/w^2+(i*k*A/(2*B))));
lag=laguerre(p,l,(k/(2*B))^2/((1/w^2)+(i*k*A/(2*B))).*(r.^2));
hgb=C1.*C2.*C3.*lag;
HGM=abs(hgb.^2);
HGM_n=HGM./max(HGM(:));
HG=HGM_n(:,500);
imagesc(X,Y,HGM_n);

colormap("hot")
shading interp; lighting phong; view(2); axis equal; axis tight;% axis off;