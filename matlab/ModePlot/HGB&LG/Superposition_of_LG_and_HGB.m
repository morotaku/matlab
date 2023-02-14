clear all
close all
clc
i=1j;
%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%%Figure parameter
width=10;
scale=1000;

%%% Topological charge
l1=2; %North pole
l2=0; %South pole
%%% Radial index
p1=0;
n=1;

%%% Poincare sphere angle
phi=0; %latitude
thita=pi/2; %longitude

%%% Polarizer angle 
angle=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%% parameter of LG mode
w1=1000; %Beam waist
lam=0.64; %Wave length
k=2*pi/lam; %Wave number
zr1=w1^2*k/2; %Rayleigh length



%%% parameter of HGB
w2=1000; %Beam waist
lam=0.63; %Wave length
zr2=w2^2*k/2; %Rayleigh length
f=100000; %Focus length
s=zr2; %Distance from the input plane to lens

w2_2=w2/(sqrt(1-(s/f)^2+(zr1/f)^2)); %Beam waist after focus
zr2_2=pi*w2_2^2/lam; %Rayleigh length after focus

z1=input('position: '); %Beam position
z2=z1+s+f;


G0=1;

A=1+(s-z2)/f;%ABCD matrix
B=z2+(s^2-2*s)/f;
C=-1/f;
D=1-s/f;

q0=1j*k*w2^2/2; % Beam q parameter
Q=(A+B/q0)/(C+D/q0);
W2=w2*sqrt(A^2+(B*lam./(pi*w2^2))^2);




%%%%%%%%%%%%%%Intensity plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1001;
L=10^4; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi1,r] = cart2pol(x,y);

%%% LG mode
LG_n=LGmode(p1,l1,r,phi1,z1,w1,lam);

%%% HGB
HGB=f_hgb(n,s,f,z2,lam,w2,w2_2,G0,r);

i_r=0;
for t = 0:20
    E1_i=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi1)).*exp(1j*t*pi/20).*LG_n;  %ℓ=1
    E2_i=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi1)).*exp(1j*t*pi/20).*HGB; %ℓ=-1
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
N2=11;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);

%%% LG mode
LG_q=LGmode(p1,l1,r2,phi2,z1,w1,lam);

%%% HGB
HGB_q=f_hgb(n,s,f,z2,lam,w2,w2_2,G0,r2);



cnt1=0;
cnt2=0;
g1=[];
g2=[];
for i=1:N2^2
    if sin(thita/2)*abs(LG_q(i))>cos(thita/2)*abs(HGB_q(i))
        cnt1=cnt1+1;
        g1(cnt1)=i;
        LG_1(cnt1)=LG_q(i);
        HGB_1(cnt1)=HGB_q(i);
        x_1(cnt1)=x2(i);
        y_1(cnt1)=y2(i);

        
    else
        cnt2=cnt2+1;
        g2(cnt2)=i;
        LG_2(cnt2)=LG_q(i);
        HGB_2(cnt2)=HGB_q(i);
        x_2(cnt2)=x2(i);
        y_2(cnt2)=y2(i);
    end
end

L2=11;
imagesc([-L L],[-L L],I);
colormap('gray')
hold on

%%%Left circular
for t = 0:120
    %%%North pole of Poincare sphere
    E1_1=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(g1))).*exp(1j*t*pi/40).*LG_1;  %ℓ=1
    %%%South pole
    E2_1=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2(g1))).*exp(1j*t*pi/40).*HGB_1; %ℓ=-1
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
for t = 0:120
    %%%North pole of Poincare sphere
    E1_2=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2(g2))).*exp(1j*t*pi/40).*LG_2;
    %%%South pole
    E2_2=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2(g2))).*exp(1j*t*pi/40).*HGB_2; %ℓ=-1
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

