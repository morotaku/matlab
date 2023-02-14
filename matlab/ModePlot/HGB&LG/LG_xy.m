clear all
close all
clc
format shortEng
%%% Topological charge
l=1;
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
LG=LGmode(p,l,r,phi,z,w,lam);
imagesc(X,Y,HGM_n);
hold on
%imagesc(X,Y,LG);



shading interp; lighting phong; view(2); axis equal; axis tight;% axis off;

