clear all
close all
clc

%%% Topological charge
l1=0; %North pole

%%% Radial index
p1=0;

%%% Polarizer angle 
angle=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%%%%%%%%%%%%%Intensity plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1001;
L=2500; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
Z=linspace(0,1.25*10^5,N);
%%% Polar coodinate
[phi1,r] = cart2pol(x,y);

%%% LG mode


for j=0:300
    z=1.25*(j/3)*10^3; %Beam position
    w=20.5; %Beam waist
    lam=0.63; %Wave length
    %z=input('position: '); %Beam position
    k=2*pi/lam; %Wave number
    zr=w^2*k/2; %Rayleigh length
    R=z+zr^2/z; %Beam curvature
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    
    lag=LGmode(p1,l1,r,phi1,z,w,lam);
    HGM=abs(lag.^2);
    HGM_n=HGM./max(HGM(:));
    %HGM_n(:,500)
    HGM_xz(:,j+1)=HGM_n(:,500);
end
    imagesc(Z,X,HGM_xz)
    colormap("hot")
    %xlim([0 4000]);
    %shading interp; lighting phong; view(2); axis equal; axis tight;
    %axis equal;
    %axis off;
