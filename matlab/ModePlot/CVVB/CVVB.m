clear all
close all
clc
%%% Topological charge
l=3;
%%% Radial index
p=0;


%%% Parameters
w=1; %Beam waist
lam=0.64; %Wave length
z=input('position: '); %Beam position
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size


%%% x-y　coordinate
N=1000;
L=10; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi,r] = cart2pol(x,y);

%%% LG mode
LG_n=LGmode(p,l,r,phi,z,w,lam);

%%% Polarizer angle 
thita=input('polarizer angle: '); %thita=-1 → no polarizer

%%% Radial mode
[ex_r,ey_r]=radial(LG_n,l,phi,thita);
i_r=realsqrt(ex_r.^2+ey_r.^2);

%%% Azimuth mode
[ex_a,ey_a]=azimuth(LG_n,l,phi,thita);
i_a=realsqrt(ex_a.^2+ey_a.^2);

%%%%%%%%%%%%%%quiver plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=20;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);
%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);
%%% LG mode
LG2_n=LGmode(p,l,r2,phi2,z,w,lam);
%%% Radial mode
[Ex_r,Ey_r]=radial(LG2_n,l,phi2,thita);
I_r=realsqrt(Ex_r.^2+Ey_r.^2);
%%% Azimuth mode
[Ex_a,Ey_a]=azimuth(LG2_n,l,phi2,thita);
I_a=realsqrt(Ex_a.^2+Ey_a.^2);


%imagesc(i_a);

imagesc([-L L],[-L L],i_r);
hold on
q=quiver(x2,y2,Ex_r,Ey_r,0.5);
q.LineWidth=1;
q.Color='red';


scale=2;
colormap("gray")
shading interp; lighting phong; view(2); axis equal; axis tight; axis off;
