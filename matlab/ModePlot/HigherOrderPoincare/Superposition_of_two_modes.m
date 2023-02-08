clear all
close all
clc

%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%%Figure parameter
width=15;
scale=0.8;

%%% Topological charge
l1=-1; %North pole
l2=0; %South pole
%%% Radial index
p1=0;
p2=0;

%%% Poincare sphere angle
phi=0; %latitude
thita=pi/2; %longitude

w=1; %Beam waist
lam=0.64; %Wave length
z=input('position: '); %Beam position
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size

%%% Polarizer angle 
angle=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%%%%%%%%%%%%%Intensity plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1000;
L=10; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi1,r] = cart2pol(x,y);

%%% LG mode
LG1_i=LGmode(p1,l1,r,phi1,z,w,lam);
LG2_i=LGmode(p2,l2,r,phi1,z,w,lam);

i_r=0;
for t = 0:20
    E1_i=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi1)).*exp(1j*t*pi/20).*LG1_i;  %ℓ=1
    E2_i=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi1)).*exp(1j*t*pi/20).*LG2_i; %ℓ=-1
    ex1_i=real(E1_i);
    ey1_i=real(-1j.*E1_i);
    ex2_i=real(E2_i);
    ey2_i=real(1j.*E2_i);
    ex_i=ex1_i+ex2_i;
    ey_i=ey1_i+ey2_i;
    [Ex_i,Ey_i]=polarizer(angle,ex_i,ey_i);
    i_r=i_r+real(sqrt(Ex_i.^2+Ey_i.^2));
end
I=i_r./max(max(i_r));

%%%%%%%%%%%%%%Quiver plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=12;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);

%%% LG mode
lg1=LGmode(p1,l1,r2,phi2,z,w,lam);
lg2=LGmode(p2,l2,r2,phi2,z,w,lam);

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


L2=11;
imagesc([-L2 L2],[-L2 L2],I);
colormap('gray')
hold on

%%%Left circular
for t = 0:60   
    %%%North pole of Poincare sphere
    E1_1=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(g1))).*exp(1j*t*pi/20).*lg1_1;  %ℓ=1
    %%%South pole
    E2_1=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2(g1))).*exp(1j*t*pi/20).*lg2_1; %ℓ=-1
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
    %pause(0.005)
end
hold on
%%%Right circular
for t = 0:60
    %%%North pole of Poincare sphere
    E1_2=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(g2))).*exp(1j*t*pi/20).*lg1_2;  %ℓ=1
    %%%South pole
    E2_2=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2(g2))).*exp(1j*t*pi/20).*lg2_2; %ℓ=-1
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
