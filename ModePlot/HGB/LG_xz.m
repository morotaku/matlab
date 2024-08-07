clear all
close all
clc
format shortEng

%%% x-y　coordinate
N=1001;
L=20; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);
Z=linspace(0,40,N);

%%% Topological charge
l=0;
%%% Radial index
p=1;


%%% Parameters(µm)
w=1; %Beam waist
lam=1.06; %Wave length
k=2*pi/lam; %Wave number

i=1j;


for j=0:1000
    z=j/10; %Beam position
    zr=w^2*k/2; %Rayleigh length
    R=z+zr^2/z; %Beam curvature
    W=w*(1+(z/zr)^2)^0.5; %Beam size 
    
    C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
    Gouy=-1j*((-(l*phi)-((2*p+abs(l)+1)*atan(z/zr)))+k*(r.^2/(2*R))+k*z);
    LG1=C/W.*exp((-1)*r.^2/W^2).*(2^0.5*r/W).^abs(l).*a_laguerre(p, l, 2*r.^2/W^2).*exp(Gouy);
    LG=abs(LG1.^2);
    LG_n=LG./max(LG(:));
    LG_z(:,j+1)=LG_n(:,500);
end
    imagesc(Z,X,LG_z)
    colormap("hot")
    shading interp; lighting phong; view(2); axis equal; axis tight; 
    %axis off;
