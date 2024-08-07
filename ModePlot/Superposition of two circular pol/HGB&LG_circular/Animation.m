clear all
close all
clc
i=1j;
L=500; %Display lange
%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%%Figure parameter
width=10;
scale=0.5;

%%% Topological charge
l1=1; %North pole
l2=0; %South pole
%%% Radial index
p1=0;
n=1;

%%% Poincare sphere angle
phi=0; %latitude
thita=pi/2; %longitude



%%% parameter of LG mode
w1=20.5; %Beam waist
lam=0.63; %Wave length
k=2*pi/lam; %Wave number
zr1=w1^2*k/2; %Rayleigh length
z1=input('position: '); %Beam position
R=z1+zr1^2/z1; %Beam curzrvature
W1=w1*(1+(z1/zr1)^2)^0.5; %Beam size

%%% parameter of HGB
w2=1000; %Beam waist
lam=0.63; %Wave length
zr2=w2^2*k/2; %Rayleigh length
f=100000; %Focus length
s=zr2; %Distance from the input plane to lens
z2=z1+f;

w2_2=w2/(sqrt(1-(s/f)^2+(zr1/f)^2));
zr2_2=pi*w2_2^2/lam;

G0=1;




%%% Polarizer angle 
angle=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%%%%%%%%%%%%%Intensity plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1001;

X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi1,r] = cart2pol(x,y);

%%% LG mode
LG_n=LGmode(p1,l1,r,phi1,z1,w1,lam);

%%% HGB
HGB=f_hgb(n,f,z2,lam,w2,G0,r);

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
HGB_q=f_hgb(n,f,z2,lam,w2,G0,r2);



for t = 0:60
    imagesc([-L L],[-L L],I);
    hold on
    %%%North pole of Poincare sphere
    E1_q=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2)).*exp(1j*t*pi/20).*LG_q;  %ℓ=1
    %%%South pole
    E2_q=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2)).*exp(1j*t*pi/20).*HGB_q; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_q=real(E1_q);
    ey1_q=real(-1j.*E1_q);
    %%%|RHC> ℓ=-1
    ex2_q=real(E2_q);
    ey2_q=real(1j.*E2_q);

    ex_q=ex1_q+ex2_q;
    ey_q=ey1_q+ey2_q;
    [Ex_q,Ey_q]=polarizer(angle,ex_q,ey_q);

    q=quiver(x2,y2,Ex_q,Ey_q,1);
    q.LineWidth=1;
    q.Color='red';
    
    colormap("gray")
    shading interp; lighting phong; view(2); axis equal; axis tight;
    hold off
    pause(0.1)
end
