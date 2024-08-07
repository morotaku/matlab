clear all
close all
clc
%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%%Figure parameter
width=10;
scale=0.75;

%%% Topological charge
l1=1; %North pole
l2=0; %South pole
%%% Radial index
p1=0;
p2=0;

L=10;
%%% Poincare sphere angle
phi=0; %latitude
thita=pi/2; %longitude


%%% Polarizer angle 
angle=0;%input('polarizer angle: '); %thita=-1 → no polarizer

gauss_pass="C:\Users\maila\OneDrive\デスクトップ\保存先\田村さん\gaussian_100";
vortex_pass="C:\Users\maila\OneDrive\デスクトップ\保存先\田村さん\vortex_100";
Gaussian=readmatrix(gauss_pass);
Vortex=readmatrix(vortex_pass);
Superposition=Gaussian+Vortex;
imagesc([-L L],[-L, L],Superposition);
colormap('gray')
axis equal; axis off;
hold on

%%%%%%%%%%%%%%Quiver plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=20;

X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);

gauss_pass="C:\Users\maila\OneDrive\デスクトップ\保存先\田村さん\gaussian_20";
vortex_pass="C:\Users\maila\OneDrive\デスクトップ\保存先\田村さん\vortex_20";
lg1=readmatrix(gauss_pass)./248;
lg2=readmatrix(vortex_pass)./248;

cnt1=0;
cnt2=0;
for i=1:N2^2
    if sin(thita/2)*abs(lg1(i))>cos(thita/2)*abs(lg2(i))
        cnt1=cnt1+1;
        g1(cnt1)=i;
        lg1_1(cnt1)=lg1(i);
        lg2_1(cnt1)=lg2(i);
        x_1(cnt1)=x2(i);
        y_1(cnt1)=y2(i);

        
    else
        cnt2=cnt2+1;
        g2(cnt2)=i;
        lg1_2(cnt2)=lg1(i);
        lg2_2(cnt2)=lg2(i);
        x_2(cnt2)=x2(i);
        y_2(cnt2)=y2(i);
    end
end

if cnt1==0
    LG_1=0;
    HGB_1=0;
    x_1=0;
    y_1=0;
    g1=1;
elseif cnt2==0
    LG_2=0;
    HGB_2=0;
    x_2=0;
    y_2=0;
    g2=1;
end

hold on


%%%Left circular
for t = 0:120   
    %%%North pole of Poincare sphere
    E1_1=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(g1))).*exp(1j*t*pi/40).*lg1_1;  %ℓ=1
    %%%South pole
    E2_1=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2(g1))).*exp(1j*t*pi/40).*lg2_1; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_1=real(E1_1);
    ey1_1=real(-1j.*E1_1);
    %%%|RHC> ℓ=-1
    ex2_1=real(E2_1);
    ey2_1=real(1j.*E2_1);
    
    ex_1=ex1_1+ex2_1;
    ey_1=ey1_1+ey2_1;
    [Ex_1,Ey_1]=polarizer(angle,ex_1,ey_1);
    Ex1=Ex_1.*scale+x_1;
    Ey1=Ey_1.*scale+y_1;
    axis equal; axis off;
    
    scatter(Ex1,Ey1,width,'.','red')
    axis equal; axis off;
    %pause(0.005)
end
hold on
%%%Right circular
for t = 0:120
    %%%North pole of Poincare sphere
    E1_2=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(g2))).*exp(1j*t*pi/40).*lg1_2;  %ℓ=1
    %%%South pole
    E2_2=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2(g2))).*exp(1j*t*pi/40).*lg2_2; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_2=real(E1_2);
    ey1_2=real(-1j.*E1_2);
    %%%|RHC> ℓ=-1
    ex2_2=real(E2_2);
    ey2_2=real(1j.*E2_2);
    
    ex_2=ex1_2+ex2_2;
    ey_2=ey1_2+ey2_2;
    [Ex_2,Ey_2]=polarizer(angle,ex_2,ey_2);
    Ex2=Ex_2.*scale+x_2;
    Ey2=Ey_2.*scale+y_2;
    axis equal; axis off;
    
    scatter(Ex2,Ey2,width,'.','blue')
    %pause(0.005)
end

function [x,y]=polarizer(thita,ex,ey)
if thita==-1
    x=ex;
    y=ey;
else
    x=cos(thita)^2.*ex+cos(thita)*sin(thita).*ey;
    y=sin(thita)*cos(thita).*ex+sin(thita)^2.*ey;
end
end

