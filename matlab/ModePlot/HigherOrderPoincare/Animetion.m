clear all
close all
clc
%%% Topological charge
l1=1; %North pole
l2=0; %South pole
%%% Radial index
p1=0;
p2=0;
%%% Poincare sphere angle
phi=0; %latitude
thita=pi/10; %longitude

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
LG1_q=LGmode(p1,l1,r2,phi2,z,w,lam);
LG2_q=LGmode(p1,l2,r2,phi2,z,w,lam);

%%% angular frequency ω
one_q=ones(size(phi2,1));
o_q=one_q.*(pi/20);


for t = 0:200
    imagesc([-L L],[-L L],i_r);
    hold on
    %%%North pole of Poincare sphere
    E1_q=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2+t*o_q)).*LG1_q;  %ℓ=1
    %%%South pole
    E2_q=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2+t*o_q)).*LG2_q; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_q=real(E1_q);
    ey1_q=real(-1j.*E1_q);
    %%%|RHC> ℓ=-1
    ex2_q=real(E2_q);
    ey2_q=real(1j.*E2_q);

    ex_q=ex1_q+ex2_q;
    ey_q=ey1_q+ey2_q;
    [Ex_q,Ey_q]=polarizer(angle,ex_q,ey_q);

    q=quiver(x2,y2,Ex_q,Ey_q,'off');
    q.LineWidth=1;
    q.Color='red';
    
    colormap("gray")
    shading interp; lighting phong; view(2); axis equal; axis tight; axis off;
    hold off
    pause(0.05)
end
