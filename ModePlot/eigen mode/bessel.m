clear all
close all
clc
%Radial index
l1=11;
l2=-l1;


%Azimuth index
p1=0;
p2=0;

%Parameters
w=1; %Beam waist
lam=0.64; %Wave length
k=2*pi/lam; %Wave number
zr=w^2*k/2; %Rayleigh length
z=input('position: '); %Beam position 
R=z+zr^2/z; %Beam curvature
W=w*(1+(z/zr)^2)^0.5; %Beam size
%Z=z/zr

%x-y　coordinate
N=1000;
L=0.1; %Display lange
X=linspace(-L,L,N);
Y=linspace(-L,L,N);
[x,y]=meshgrid(X,Y);
%Polar coodinate
%r=sqrt(x^2+y^2);
%phi=atan2(y,x);
[phi,r] = cart2pol(x,y);
ang=pi/18;
l1=7;
l2=-7;
kr=sqrt(tan(ang)^2/(1+tan(ang)^2))*(pi/lam^2);
kz=sqrt(4*pi/lam^2-kr^2);
a=0.1;
alf=570./(12.*a.*z).*(1+10.*a.*z+sqrt(1+20.*a.*z+4.*a.^2.*z.^2));
BBeam1=besselj(l1,kr.*r.*alf).*exp(1j.*l1.*phi).*exp(1j.*kz.*z);
BBeam2=besselj(l2,kr.*r.*alf).*exp(1j.*l2.*phi).*exp(1j.*kz.*z);
BBeam=BBeam1+BBeam2;
%LG=abs((LG1+LG2).^2);
BesselBeam=abs((BBeam).*conj(BBeam));
BesselBeam_n=BesselBeam./max(BesselBeam(:));
BesselBeam_a=angle(BBeam);
imagesc(BesselBeam_n,[0 1]);

G=[1 1 1];
for i=0:50
    G=[G;[1-0.02*i 1 1-0.02*i]];
end
colormap(gray)
shading interp; lighting phong; view(2); axis equal; axis tight; axis off;

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
%function y=LGMode(p,l,rad)
%C=(2*factorial(p)./(pi*factorial(p+abs(l))))^0.5;
%Gouy=-1j*((-(l.*phi)-((2*p+abs(l)+1)*atan(z/zr)))+k.*(rad.^2./(2*R))+k*z);
%y=C/W.*(2^0.5.*rad./W)^abs(l).*laguerre(p, l, 2.*rad^2./W^2).*exp((-1).*rad.^2./W^2).*exp(Gouy);
%end
