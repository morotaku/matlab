clear all
close all
clc
i=1j;
L=9;
pass='C:\Users\maila\OneDrive\デスクトップ\保存先\4\270.emf';
%%%Figure parameter
width=15;
scale=L/3;
ang=input('angle=');
%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%

%%% Topological charge
l1=0; %North pole
l2=0; %South pole
%%% Radial index
p1=0;
p2=0;

%%% Poincare sphere angle
phi=pi; %latitude
thita=3*pi/4; %longitude強度比

%%% Polarizer angle 
angle=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%% parameter of LG mode
w1=5; %Beam waist
lam=0.63; %Wave length
k=2*pi/lam; %Wave number
zr1=w1^2*k/2 %Rayleigh length



%%% parameter of HGB
w2=5; %Incident Beam size
lam=0.63; %Wave length
zr2=w2^2*k/2; %Rayleigh length
f=10^5; %Focus length
s=zr2; %Distance from the input plane to lens

w2_2=w2/(sqrt(1-(s/f)^2+(zr1/f)^2)); %Beam waist after focus
zr2_2=pi*w2_2^2/lam; %Rayleigh length after focus

z1=input('position: '); %Beam position
z2=z1;


G0=1;


%%%%%%%%%%%%%%Intensity plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N=1001;
%L=35; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);

%%% Polar coodinate
[phi1,r] = cart2pol(x,y);

%%% LG mode
LG_i1=LGmode(p1,l1,r,phi1,z1,w1,lam);

%%% HGB
LG_i2=LGmode(p2,l2,r,phi1,z1,w2,lam);

i_x=0; i_y=0;
i_r=0;
for t = 0:40
    E1_i=sin(thita/2)*exp(1j*phi/2).*exp(1j*t*pi/20).*LG_i1;  %ℓ=1
    E2_i=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/20).*LG_i2.*0; %ℓ=-1
    ex1_i=E1_i;
    ey1_i=-1j.*E1_i.*0;
    ex2_i=E2_i;
    ey2_i=i.*E2_i;
    ex_i=ex1_i+ex2_i;
    ey_i=ey1_i+ey2_i;
    [Ex_i,Ey_i]=polarizer(pi/2,ex_i,ey_i);
    %i_r=i_r+real((Ex_i.*conj(Ex_i)+Ey_i.*conj(Ey_i)));
    i_x=i_x+Ex_i;
    i_y=i_y+Ey_i;
end
I=real(i_x.*conj(i_x))+real(i_y.*conj(i_y));
I_n=I./max(max(I));
%I=i_r./max(max(i_r));

%axis off;
figure(1)
imagesc([-L L],[-L L],I_n);
xticks(-L:L/5:L)
yticks(-L:L/5:L)
colormap('gray')
axis equal; axis off;
shading interp; lighting phong; view(2); axis tight;
hold on
%%%%%%%%%%%%%%Quiver plot%%%%%%%%%%%%%%%%%%%%
%%% x-y　coordinate
N2=15;
X2=linspace(-L,L,N2);
Y2=linspace(-L,L,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);

%%% LG mode
LG_r=LGmode(p1,l1,r2,phi2,z1,w1,lam);

%%% HGB
LG_l=LGmode(p1,l2,r2,phi2,z1,w2,lam);


%%%Right circular
for t = 0:120
    %%%North pole of Poincare sphere
    Er=sin(thita/2)*exp(1j*phi/2).*exp(1j*t*pi/60).*LG_r;  %ℓ=1
    %%%South pole
    El=cos(thita/2)*exp(-1j*phi/2).*exp(1j*t*pi/60).*LG_l; %ℓ=-1
    %%%|LHC> ℓ=1
    ex_r=real(Er);
    ey_r=real(-1j.*Er);
    %%%|RHC> ℓ=-1
    ex_l=real(El);
    ey_l=real(1j.*El);
    
    ex=ex_r+ex_l;
    ey=ey_r+ey_l;
    [Ex,Ey]=polarizer(angle,ex,ey);
    Ex_i=Ex.*scale+x2;
    Ey_i=Ey.*scale+y2;
    axis equal; axis off;
    
    scatter(Ex_i,Ey_i,width,'.','green')
    pause(0.001)
end



function y=LGmode(p,l,r,phi,z,w,lam)
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2./z; %Beam curvature
W=w.*(1+(z./zr).^2).^0.5; %Beam size
C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
Gouy=-1j.*((-(l*phi)-((2*p+abs(l)+1).*atan(z./zr)))+k.*(r.^2./(2*R))+k.*z);
lg=C./W.*exp((-1).*r.^2./W.^2).*(2^0.5*r./W).^abs(l).*laguerre(p, l, 2*r.^2./W.^2).*exp(Gouy);
LG=abs(lg.^2);
%y=LG./max(LG(:));
y=lg;
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

function [x,y]=polarizer(thita,ex,ey)
if thita==-1
    x=ex;
    y=ey;
else
    x=cos(thita)^2.*ex+cos(thita)*sin(thita).*ex;
    y=sin(thita)*cos(thita).*ex+sin(thita)^2.*ey;
end
end

