clear all
close all
clc

color='jet';
%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%%%%%%%%%
%%% Topological charge
l1=1; 
l2=-1; 
%%% Radial index
p1=0;
p2=0;

%%% Input beam intensity E1:LG(1,0) E2:LG(-1,0)
Ex1=1;
Ex2=-Ex1; %Ex2=Ex1 → Ex=HG01 , Ex2=-Ex1 → Ex=HG10 
Ey1=0; %Ey1=1j→phase defference pi/2 between x- and y-polarization 
Ey2=-Ey1; %Ey2=Ey1 → Ey=HG01 , Ey2=-Ey1 → Ey=HG10

%%% Polarizer angle 
ang=-1;%input('polarizer angle: '); %thita=-1 → no polarizer

%%% parameter of input beam
w1=1000; %Beam waist
w2=w1;
lam=1.064; %Wave length
k=2*pi/lam; %Wave number
%%% Parameter
NA=0.96;
f=1.7*10^3;
w0=0.61*lam/NA
zr=pi*w0^2/lam
z=zr.*0;
alpha=asin(NA);


%z=input('position: '); %Beam position

%% Input beam image
%%% x-y　coordinate
N=301;
L=2;
width=10;
scale=0.005;

X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
%%% Polar coodinate
[phi,r] = cart2pol(x,y);


%% Intensity plot
%%%Debye-Wolf integral
eX_x1=Ex1.*1j^l1.*(I0(l1,l1,alpha,k,z,phi,r,f)+I1(l1,l1+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)+I1(l1,l1-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eX_x2=Ex2.*1j^l2.*(I0(l2,l2,alpha,k,z,phi,r,f)+I1(l2,l2+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)+I1(l2,l2-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eX_y1=Ey1.*1j^(l1).*(-I1(l1,l1+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)+I1(l1,l1-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eX_y2=Ey2.*1j^(l2).*(-I1(l2,l2+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)+I1(l2,l2-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eX=(eX_x1+eX_x2)+(eX_y1+eX_y2);

eY_x1=Ex1.*1j^(l1+1).*(-I1(l1,l1+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)+I1(l1,l1-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eY_x2=Ex2.*1j^(l2+1).*(-I1(l2,l2+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)+I1(l2,l2-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eY_y1=Ey1.*1j^(l1-1).*(I0(l1,l1,alpha,k,z,phi,r,f)-I1(l1,l1+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)-I1(l1,l1-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eY_y2=Ey2.*1j^(l2-1).*(I0(l2,l2,alpha,k,z,phi,r,f)-I1(l2,l2+2,alpha,k,z,phi,r,f).*exp(2.*1j.*phi)-I1(l2,l2-2,alpha,k,z,phi,r,f).*exp(-2.*1j.*phi));
eY=(eY_x1+eY_x2)+(eY_y1+eY_y2);

eZ_x1=Ex1.*1j.^(l1-1).*(I2(l1,l1+1,alpha,k,z,phi,r,f).*exp(1j.*phi)-I2(l1,l1-1,alpha,k,z,phi,r,f).*exp(-1j.*phi));
eZ_x2=Ex2.*1j.^(l2-1).*(I2(l2,l2+1,alpha,k,z,phi,r,f).*exp(1j.*phi)-I2(l2,l2-1,alpha,k,z,phi,r,f).*exp(-1j.*phi));
eZ_y1=Ey1.*1j.^(l1-1).*(-I2(l1,l1+1,alpha,k,z,phi,r,f).*exp(1j.*phi)-I2(l1,l1-1,alpha,k,z,phi,r,f).*exp(-1j.*phi));
eZ_y2=Ey2.*1j.^(l2-1).*(-I2(l2,l2+1,alpha,k,z,phi,r,f).*exp(1j.*phi)-I2(l2,l2-1,alpha,k,z,phi,r,f).*exp(-1j.*phi));
eZ=(eZ_x1+eZ_x2)+(eZ_y1+eZ_y2);


iX=eX.*conj(eX);
iY=eY.*conj(eY);
iZ=eZ.*conj(eZ);
it=iX+iY+iZ;

IX=iX./max(max(it));
IY=iY./max(max(it));
IZ=iZ./max(max(it));
IT=it./max(max(it));


f1=figure(1);
imagesc(X,Y,IT);
title('total intensity')
xticks(-L:L/2:L);yticks(-L:L/2:L)
xlim([-L L]);ylim([-L L])
fontsize(1,25,"points")
axis tight;axis equal;
colormap('gray')
colorbar
hold on


%% polarization distribution plot
%%% x-y　coordinate
N2=15;
L2=2;
X2=linspace(-L2,L2,N2);
Y2=linspace(-L2,L2,N2);
[x2,y2]=meshgrid(X2,Y2);

%%% Polar coodinate
[phi2,r2] = cart2pol(x2,y2);


eX_x1=Ex1.*1j^l1.*(I0(l1,l1,alpha,k,z,phi2,r2,f)+I1(l1,l1+2,alpha,k,z,phi2,r2,f).*exp(2.*1j.*phi2)+I1(l1,l1-2,alpha,k,z,phi2,r2,f).*exp(-2.*1j.*phi2));
eX_x2=Ex2.*1j^l2.*(I0(l2,l2,alpha,k,z,phi2,r2,f)+I1(l2,l2+2,alpha,k,z,phi2,r2,f).*exp(2.*1j.*phi2)+I1(l2,l2-2,alpha,k,z,phi2,r2,f).*exp(-2.*1j.*phi2));
eX_y1=Ey1.*1j^l1.*(-I1(l1,l1+2,alpha,k,z,phi2,r2,f).*exp(2.*1j.*phi2)+I1(l1,l1-2,alpha,k,z,phi2,r2,f).*exp(-2.*1j.*phi2));
eX_y2=Ey2.*1j^l2.*(-I1(l2,l2+2,alpha,k,z,phi2,r2,f).*exp(2.*1j.*phi2)+I1(l2,l2-2,alpha,k,z,phi2,r2,f).*exp(-2.*1j.*phi2));
eX=(eX_x1+eX_x2)+(eX_y1+eX_y2);

eY_x1=Ex1.*1j^(l1+1).*(-I1(l1,l1+2,alpha,k,z,phi2,r2,f).*exp(2.*1j.*phi2)+I1(l1,l1-2,alpha,k,z,phi2,r2,f).*exp(-2.*1j.*phi2));
eY_x2=Ex2.*1j^(l2+1).*(-I1(l2,l2+2,alpha,k,z,phi2,r2,f).*exp(2.*1j.*phi2)+I1(l2,l2-2,alpha,k,z,phi2,r2,f).*exp(-2.*1j.*phi2));
eY_y1=Ey1.*1j^(l1-1).*(I0(l1,l1,alpha,k,z,phi2,r2,f)-I1(l1,l1+2,alpha,k,z,phi2,r2,f).*exp(2.*1j.*phi2)-I1(l1,l1-2,alpha,k,z,phi2,r2,f).*exp(-2.*1j.*phi2));
eY_y2=Ey2.*1j^(l2-1).*(I0(l2,l2,alpha,k,z,phi2,r2,f)-I1(l2,l2+2,alpha,k,z,phi2,r2,f).*exp(2.*1j.*phi2)-I1(l2,l2-2,alpha,k,z,phi2,r2,f).*exp(-2.*1j.*phi2));
eY=(eY_x1+eY_x2)+(eY_y1+eY_y2);

eZ_x1=Ex1.*1j.^(l1-1).*(I2(l1,l1+1,alpha,k,z,phi2,r2,f).*exp(1j.*phi2)-I2(l1,l1-1,alpha,k,z,phi2,r2,f).*exp(-1j.*phi2));
eZ_x2=Ex2.*1j.^(l2-1).*(I2(l2,l2+1,alpha,k,z,phi2,r2,f).*exp(1j.*phi2)-I2(l2,l2-1,alpha,k,z,phi2,r2,f).*exp(-1j.*phi2));
eZ_y1=Ey1.*1j.^(l1-1).*(-I2(l1,l1+1,alpha,k,z,phi2,r2,f).*exp(1j.*phi2)-I2(l1,l1-1,alpha,k,z,phi2,r2,f).*exp(-1j.*phi2));
eZ_y2=Ey2.*1j.^(l2-1).*(-I2(l2,l2+1,alpha,k,z,phi2,r2,f).*exp(1j.*phi2)-I2(l2,l2-1,alpha,k,z,phi2,r2,f).*exp(-1j.*phi2));
eZ=(eZ_x1+eZ_x2)+(eZ_y1+eZ_y2);

for t = 0:120
    ex=real(exp(1j*t*pi/60).*eX);  %ℓ=1
    ey=real(exp(1j*t*pi/60).*eY); %ℓ=-1

    [Ex_1,Ey_1]=polarizer(ang,ex,ey);
    Ex1=Ex_1.*scale+x2;
    Ey1=Ey_1.*scale+y2;
    axis equal; %axis off;
    
    scatter(Ex1,Ey1,width,'.','red')
    pause(0.005)
end
shading interp; lighting phong; view(2); axis tight;
axis off;









%% Function
%%%Debye-Wolf integral part1
function y=I0(l1,l2,alpha,k,z,phi,r,f)
    I=0;
    w0=1000;lam=1.064;
    for t=0:30
        w=alpha/30;
        thita=w*t;
        r_int=f.*tan(thita);
        A=LGmode(0,l1,r_int,phi,z,w0,lam);
        I=I+sqrt(cos(thita)).*sin(thita).*(1+cos(thita)).*besselj(l2,k.*r.*sin(thita)).*exp(1j.*k.*z.*cos(thita)).*A;
    end
    y=-1j*k.*f.*I;
end
%%%Debye-Wolf inegral part2
function y=I1(l1,l2,alpha,k,z,phi,r,f)
    I=0;
    w0=1000;lam=1.064;
    for t=0:30
        w=alpha/30;
        thita=w*t;
        r_int=f.*tan(thita);
        A=LGmode(0,l1,r_int,phi,z,w0,lam);
        I=I+sqrt(cos(thita)).*sin(thita).*(1-cos(thita)).*besselj(l2,k.*r.*sin(thita)).*exp(1j.*k.*z.*cos(thita)).*A;
    end
    y=-1j*k.*f.*I/2;
end
%%%Debye-Wolf integral part3
function y=I2(l1,l2,alpha,k,z,phi,r,f)
    I=0;
    w0=1000;lam=1.064;
    for t=0:30
        w=alpha/30;
        thita=w*t;
        r_int=f.*tan(thita);
        A=LGmode(0,l1,r_int,phi,z,w0,lam);
        I=I+sqrt(cos(thita)).*sin(thita).^2.*besselj(l2,k.*r.*sin(thita)).*exp(1j.*k.*z.*cos(thita)).*A;
    end
    y=-1j*k.*f.*I;
end
%%%LG
function y=LGmode(p,l,r,phi,z,w,lam)
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
R=z+zr^2./z; %Beam curvature
W=w.*(1+(z./zr).^2).^0.5; %Beam size
C=(2*factorial(p)/(pi*factorial(p+abs(l))))^0.5;
%Gouy=-1j.*((-(l*phi)-((2*p+abs(l)+1).*atan(z./zr)))+k.*(r.^2./(2*R))+k.*z);
Gouy=-1j.*((-(l*phi)-((2*p+abs(l)+1).*atan(z./zr)))+k.*z);
LG1=C./W.*exp((-1).*r.^2./W.^2).*(2^0.5*r./W).^abs(l).*laguerre(p, l, 2*r.^2./W.^2).*exp(Gouy);
LG=abs(LG1.^2);
%y=LG./max(LG(:));
y=LG1;
end
%%%associated laguerre polynomial
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
%%% HG
function y=HGmode(n,m,x,y,w,r,z,R,k)
    zr=w^2*k/2; %Rayleigh length
    gouy=1j*(m+n+1)*atan(z/zr);
    W=w*(1+(z/zr)^2)^0.5; %Beam size
    u=(1/W).*HermitePol(n,sqrt(2)/W.*x).*HermitePol(m,sqrt(2)/W.*y).*exp(-(r.^2./W^2)-(1j*k.*(x.^2+y.^2)./R));%.*exp(gouy);
    U=abs(u.^2);
    U_n=U./max(U(:));
    y=u;%./max(max(abs(u)));
end
%%% Hermite polynomial
function y=HermitePol(n,x)
Hn=0;
for m=0:floor(n/2)
    Hn=Hn+((-1)^m.*(2.*x).^(n-2*m))./(factorial(m)*factorial(n-2*m));
end
Hn=Hn.*factorial(n);
y=Hn;
end
%%% polarizert
function [x,y]=polarizer(thita,ex,ey)
if thita==-1
    x=ex;
    y=ey;
else
    x=cos(thita)^2.*ex+cos(thita)*sin(thita).*ey;
    y=sin(thita)*cos(thita).*ex+sin(thita)^2.*ey;
end
end