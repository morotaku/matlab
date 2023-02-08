clear all
close all
clc
map=[0 1 0];
%%% Topological charge
l1=1; %North pole
l2=0; %South pole
%%% Radial index
p1=0;
p2=0;
%%% Poincare sphere angle
phi=3*pi/2; %latitude
thita=pi/8; %longitude

%%% Parameters
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
E1_i=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi1)).*LG1_i; % ℓ=ℓ
E2_i=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi1)).*LG2_i; % ℓ=-ℓ

%%%Linear → Circular
ex1_i=real(E1_i); 
ey1_i=real(-1j.*E1_i);
ex2_i=real(E2_i);
ey2_i=real(1j.*E2_i);
ex_i=ex1_i+ex2_i;
ey_i=ey1_i+ey2_i;
[Ex_i,Ey_i]=polarizer(angle,ex_i,ey_i);
i_r=real(sqrt(Ex_i.^2+Ey_i.^2));

%%%%%%%%%%%%%%Quiver plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=10;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);

%%% LG mode
lg1=LGmode(p1,l1,r2,phi2,z,w,lam);
lg2=LGmode(p2,l2,r2,phi2,z,w,lam);

%lg1_1=zeros(N2,N2);
%lg2_1=zeros(N2,N2);
%lg1_2=zeros(N2,N2);
%lg2_2=zeros(N2,N2);
cnt1=0;
cnt2=0;
for i=1:100
    if lg1(i)>=lg2(i)
        cnt1=cnt1+1;
        lg1(cnt1)=i;
        lg1_1(cnt1)=lg1(i);
        lg2_1(cnt1)=lg2(i);
        x_1(cnt1)=x2(i);
        y_1(cnt1)=y2(i);
        %lg1_2(i)=0;
        %lg2_2(i)=0;

        
    else
        cnt2=cnt2+1;
        lg2(cnt2)=i;
        lg1_2(cnt2)=lg1(i);
        lg2_2(cnt2)=lg2(i);
        x_2(cnt2)=x2(i);
        y_2(cnt2)=y2(i);
        %lg1_2(i)=0;
        %lg1_1(i)=0;
        %lg2_1(i)=0;
        %cnt2=cnt2+1;
        %lg1(cnt2)=i;
    end
end
