clear all
close all
clc

%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%%Figure parameter
width=10;
scale=0.75;

%%% Topological charge
l1=1; %North pole
l2=-1; %South pole
%%% Radial index
p1=0;
p2=0;
%%% Polarizer angle 
angle=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%% Poincare sphere angle
phi=0; %latitude
thita=pi/2; %longitude

w=1; %Beam waist
lam=0.64; %Wave length
z=input('position: '); %Beam position
k=2*pi/lam; %Wave number
zr=w^2*k/2 %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size


%%%%%%%%%%%%%%Intensity plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1000;
L=8; %Display lange
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
    ey1_i=0;
    ex2_i=0;
    ey2_i=real(E2_i);
    ex_i=ex1_i+ex2_i;
    ey_i=ey1_i+ey2_i;
    [Ex_i,Ey_i]=polarizer(angle,ex_i,ey_i);
    i_r=i_r+real(sqrt(Ex_i.^2+Ey_i.^2));
end
I=i_r./max(max(i_r));
imagesc([-L L],[-L L],I);
colormap('gray')
hold on
axis equal; axis off;
%%%%%%%%%%%%%%Quiver plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=11;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);

%%% LG mode
lg1=LGmode(p1,l1,r2,phi2,z,w,lam);
lg2=LGmode(p2,l2,r2,phi2,z,w,lam);


%%%Linear
for t = 0:120
    %%%North pole of Poincare sphere
    E1_1=sin(thita/2)*exp(1j*phi/2).*exp(1j.*(l1.*phi2)).*exp(1j*t*pi/60).*lg1;  %ℓ=1
    %%%South pole
    E2_1=cos(thita/2)*exp(-1j*phi/2).*exp(1j.*(l2.*phi2)).*exp(1j*t*pi/60).*lg2; %ℓ=-1
    %%%|LHC> ℓ=1
    ex1_1=real(E1_1);
    ey1_1=0;
    %%%|RHC> ℓ=-1
    ex2_1=0;
    ey2_1=real(E2_1);
    
    ex_1=ex1_1+ex2_1;
    ey_1=ey1_1+ey2_1;
    [Ex_1,Ey_1]=polarizer(angle,ex_1,ey_1);
    Ex1=Ex_1.*scale+x2;
    Ey1=Ey_1.*scale+y2;
    axis equal; %axis off;
    
    scatter(Ex1,Ey1,width,'.','red')
    pause(0.005)
end
shading interp; lighting phong; view(2); axis tight;
axis off;


function y=LGmode(p,l,r,phi,z,w,lam)
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size
C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
Gouy=-1j*((-(l*phi)-((2*p+abs(l)+1)*atan(z/zr)))+k*(r.^2/(2*R))+k*z);
LG1=C/W.*exp((-1)*r.^2/W^2).*(2^0.5*r/W).^abs(l).*laguerre(p, l, 2*r.^2/W^2).*exp(Gouy);
LG=abs(LG1.^2);
LG_n=LG./max(LG(:));
y=LG_n;
end

function y=laguerre(p,l,x)
y=zeros(p+1,1);
L=abs(l);
if p==0
    y=1;
else
for m=0:p
    y(p+1-m)=((-1).^m.*(factorial(p+L)))./(factorial(p-m).*factorial(L+m).*factorial(m));
end
end
y=polyval(y,x);
end
