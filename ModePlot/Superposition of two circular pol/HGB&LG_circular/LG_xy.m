clear all
close all
clc
format shortEng
%%% Topological charge
l=0;
%%% Radial index
p=0;

%%% Parameters(µm)
w=1000; %Beam waist
lam=0.63; %Wave length
k=2*pi/lam; %Wave number
zr=w^2*k/2 %Rayleigh length
z=input('position: ');


%%%%%%%%%%%%%%%%Intinsity plot%%%%%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1001;
L=10000; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
[phi,r] = cart2pol(x,y);

%%%Mode plot
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size
C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
Gouy=-1j*((-(l*phi)-((2*p+abs(l)+1)*atan(z/zr)))+k*(r.^2/(2*R))+k*z);
LG1=C/W.*exp((-1)*r.^2/W^2).*(2^0.5*r/W).^abs(l).*laguerre(p,2*r.^2/W^2).*exp(Gouy);
LG=abs(LG1.^2);
LG_n=LG./max(LG(:));
imagesc(X,Y,LG_n);
hold on
%imagesc(X,Y,LG);
%colormap("hot")
Red=[1 1 1];
for i=0:50
    Red=[Red;[1 1-0.02*i 1-0.02*i]];
end
%mycolors = [1 1 1; 1 . 0];
colormap(Red)


shading interp; lighting phong; view(2); axis equal; axis tight; axis off;

