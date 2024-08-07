clear all
close all
clc
%%% Topological charge
l1=0; %topological charge of HGB
l2=1; %topological charge of LG
%%% Radial index
p1=1; %Radial index of HGB
p2=0; %Radial index of LG
%%% Poincare sphere angle
phi=0;
thita=pi/2;

z=input('position: '); %Beam position
lam=1.064; %Wave length
k=2*pi/lam; %Wave number

%%% Parameters of HGB
w1=10000; %Beam waist at incident plane
zr1=w1^2*k/2; %Rayleigh length
R1=z+zr1^2/z; %Beam curvature
W1=w1*(z+(z/zr1)^2)^0.5; %Beam size

%%% Parameters of LG
w2=8.4715; %Beam waist
zr2=w2^2*k/2; %Rayleigh length
R2=z+zr2^2/z; %Beam curvature
W2=w2*(z+(z/zr2)^2)^0.5; %Beam size



f=20000; %Focus length
s=200000; %Distance from the input plane to lens
G0=1;
i=1j;

%ABCD matrix
A=-z/f;
B=-(z/f)*s+f+z;
C=-1/f;
D=1-s/f;

%%% Polarizer angle 
angle=input('polarizer angle: '); %thita=-1 → no polarizer


%%%%%%%%%%%%%%Intensity plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N1=1001;
L=20; %Display lange
L=20;
X1=linspace(-L,L,N1);
Y1=linspace(-L,L,N1);
[x1,y1]=meshgrid(X1,Y1);
[phi1,r1] = cart2pol(x1,y1);

%%% LG mode
HGB_i=hgb(p1,l1,r1,z,w1,k,A,B,C,D,G0);
LG_i=LGmode(p2,l2,r1,phi1,0,w2,W2,lam,R2);
E1_i=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi1)).*HGB_i; % ℓ=ℓ
E2_i=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(-l2.*phi1)).*LG_i; % ℓ=-ℓ
HGB_i(:,500);
e=2.71;

for i=500:1001
    if HGB_i(1001-i,500)>(1/2)
        FWHD=(L/N1)*(1001-i)
        break
    end
end

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
N2=20;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);

%%% LG mode
HGB_q=hgb(p1,l1,r2,z,w1,k,A,B,C,D,G0);
LG_q=LGmode(p2,l2,r2,phi2,0,w2,W2,lam,R2);

%%% angular frequency ω
one_q=ones(size(phi2,1));
o_q=one_q.*(pi/20);

for t = 0:200
    imagesc(X1,Y1,HGB_i);
    hold on
    %%%North pole of Poincare sphere
    E1_q=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2+t*o_q)).*HGB_q;  %ℓ=1
    %%%South pole
    E2_q=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(-l2.*phi2+t*o_q)).*LG_q; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_q=real(E1_q);
    ey1_q=real(-1j.*E1_q);
    %%%|RHC> ℓ=-1
    ex2_q=real(E2_q);
    ey2_q=real(1j.*E2_q);

    ex_q=ex1_q+ex2_q;
    ey_q=ey1_q+ey2_q;
    [Ex_q,Ey_q]=polarizer(angle,ex_q,ey_q);

    q=quiver(x2,y2,ex1_q,ey1_q);
    q.LineWidth=1;
    q.Color='red';
    colormap("gray")
    shading interp; lighting phong; view(2); axis equal; axis tight;% axis off;
    hold off
    pause(0.05)
end
